#!/usr/bin/perl
# Wed 21 Mar 2012 02:30:34 PM PDT automates initial library prep
# Fri 23 Mar 2012 10:17:03 AM PDT modified to work with a pair of libraries
# Wed 04 Jul 2012 08:07:18 AM PDT modified to work on hbar, removed Temp:seekable dependance and prevented making graphs because system lacks "gnuplot" library
# Fri 06 Jul 2012 10:57:13 AM PDT added a test to see if data is ok for chastity filter
# Sept/Oct 2012 modified the ribosome matching to accomodate spaces
# Fri 05 Oct 2012 12:59:25 PM PDT changed output to comma separated
# April 2014 now using SortMeRNA as ribosome removal program, updated Trimmomatic
# Adapted for cleanning database prior to RepeatExplorer assembly

use strict;
use File::Temp ();
use File::Basename;
use Getopt::Long;
use File::Path;

# CONSTANTS
# adapter trimming parameters
my $TRIMMOMATIC_PATH = "/home/parensburge/RNAseq-clean/Trimmomatic-0.32"; # java program location
# my $SORTMELOC = "/home/parensburge/Desktop/RNAseq-clean/sortmerna-1.99-beta-linux-64bit"; # must keep full path name for sortme program
my $MINLEN = 35;

#return date and time
sub datetime {
	use POSIX qw/strftime/;
	return (strftime('%D %T',localtime));
}

# other parameters
my $readspair1; #filename of first pair fastq
my $readspair2; #filename of second pair fastq file
my $outputname; #base name for output files
my $outputdir; #outputdirectory
my $threads = `grep -c processor /proc/cpuinfo`; #number of threads to use
$threads =~ s/\s//g;

#####read and check the inputs
GetOptions(
	'1:s'   => \$readspair1,
	'2:s'	=> \$readspair2,
	'o:s'	=> \$outputname,
	'd:s'	=> \$outputdir,
	't:s'   => \$threads,
);

unless ($readspair1 and $readspair2) {
	die ("usage: perl clean_pair.pl -1 <FASTQ file pair 1> -2 <FASTQ file pair 2> -o <OPTIONAL output name> -d <OPIONAL output directory (default current directory)> -t <OPTIONAL: number of threads to use, default $threads>\n");
}
unless ($outputname) {
	$outputname = basename($readspair1, ".fastq"); #take the name from first member of the pair
}

#test if output directory has been specified, if not set output to current directoy
unless ($outputdir) {
	$outputdir = `pwd`; #set output to current directoy
	chomp $outputdir;
}

#test if output directory exists if not, create it
unless ( -e $outputdir )
{
	my $dir_test = mkdir($outputdir, 0777);
	unless ($dir_test) { 
		die "cannot create directory $outputdir\n";
	}
}

#create temporary files
my $paired_output = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # temporary file with results pair1
my $unpaired_output = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # temporary file with results pair1

#create log file
if ( -e $outputdir) {
	open (LOG, ">$outputdir/Log.txt") or die ("cannot create $outputdir/Log.txt");
}
else {
	open (LOG, ">Log.txt") or die ("cannot create Log.txt");
}

########### start cleanning ######################
print LOG datetime, " Initial count\n";
print LOG datetime, " File $readspair1, total FASTQ reads: ", count_fastq($readspair1), "\n"; # do a basic count
print LOG datetime, " File $readspair2, total FASTQ reads: ", count_fastq($readspair2), "\n"; # do a basic count
print LOG "\n";

#clip adapters
print LOG datetime, " Clipping adapters and quality filter with Trimmomatic\n";
clipadapters($readspair1, $readspair2);
print LOG datetime, " File with paired reads, FASTQ reads: ", count_fastq($paired_output), "\n"; 
print LOG datetime, " File with unpaired reads, FASTQ reads: ", count_fastq($unpaired_output), "\n";
print LOG "\n";

# filter by quality
my $min_qual = 10;
my $percent = 90; # ratio of high quality sequences over length of sequence
print LOG datetime, " Filtering low quality reads using fastq_quality_filter with -q $min_qual -p $percent\n";
filterbyquality($paired_output, $min_qual, $percent);
print LOG datetime, " File with paired reads, FASTQ reads: ", count_fastq($paired_output), "\n";

# filter removing overlapping reads
my $min_overlap = 30;
my $max_overlap = 250;
print LOG datetime, " Removing read pairs that overlap by at least $min_overlap bp.\n";
remove_overlaps($paired_output, $min_overlap, $max_overlap);
print LOG datetime, " File with paired reads, FASTQ reads: ", count_fastq($paired_output), "\n";

#print the data files
my $pairedoutname = $outputdir . "/" . $outputname .  "-paired.fq"; 
`mv $paired_output $pairedoutname`; 
print LOG datetime, " Paired data file is written in $pairedoutname\n";
print LOG "\n";

######## subroutines ###################
sub count_fastq {
#	print STDERR "count fastq...\n";
	my ($title) = @_;
	my $rawcount;

	my $txtcount = `wc -l $title`;
	if ($txtcount =~ /^(\d+)\s/) {
		my $count = $1/4;
		return (commify($count));
#		return(sprintf('%.3e', $count));
#		return ($count);
	}
	else {
		die "could not count lines in file $title using wc -l\n";
	}
}

sub filterbyquality {
	my ($inputfile, $min_qual, $percent) = @_;
	my $tempout = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # file that has unpaired sequences
	
	`fastq_quality_filter -q $min_qual -p $percent -i $inputfile -o $tempout`; 
	# remove unpaired sequences and write out paired sequences
	open (OUTPUT, ">$paired_output") or die;
	open (INPUT, $tempout) or die "cannot open file $tempout\n";
	(my $l1) = split /[\s]+/, <INPUT>, 2; # this matches only the first word of the string, got this from http://stackoverflow.com/questions/4973229/perl-get-first-word-from-input-string
	my $l2 = <INPUT> . <INPUT> . <INPUT>;
	while (my $l3 = <INPUT>) {
		($l3) = split /[\s]+/,  $l3, 2;
		my $l4 = <INPUT> . <INPUT> . <INPUT>;
		if ($l1 eq $l3) {
			print OUTPUT "$l1\n", $l2, "$l3\n", $l4;
			($l1) = split /[\s]+/, <INPUT>, 2;
			$l2 = <INPUT> . <INPUT> . <INPUT>;
		}
		else {
			$l1 = $l3;
			$l2 = $l4;
		}	
	}
	close INPUT;
	close OUTPUT;
}

sub clipadapters {
#	print STDERRR "clipping adapters...\n";	
	my ($inputfile1, $inputfile2) = @_;
	my $forward_paired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); 
	my $reverse_paired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); 
	my $forward_unpaired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' );
	my $reverse_unpaired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); 

	`java -jar $TRIMMOMATIC_PATH/trimmomatic-0.32.jar PE -threads $threads -phred33 $inputfile1 $inputfile2 $forward_paired $forward_unpaired $reverse_paired $reverse_unpaired ILLUMINACLIP:$TRIMMOMATIC_PATH/illuminaClipping.fa:2:30:10 `;
#SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:$MINLEN	
	#merge the files
	open (INPUT1, $forward_paired) or die;
	open (INPUT2, $reverse_paired) or die;
	open (OUTPUT1, ">$paired_output") or die;
	while (my $line1 = <INPUT1>) {
		my $seq1 = $line1 . <INPUT1> . <INPUT1> . <INPUT1>;
		my $seq2 = <INPUT2> . <INPUT2> . <INPUT2> . <INPUT2>;
		print OUTPUT1 "$seq1";
		print OUTPUT1 "$seq2";
	}	
	`cat $forward_unpaired $reverse_unpaired > $unpaired_output`;
	close INPUT1;
	close INPUT2;
	close OUTPUT1; 
}

# from perl cookbook, to put commas on big numbers
sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text
}

sub remove_overlaps {
	my ($input_file, $min_overlap, $max_overlap) = @_;
	my $tempdir = File::Temp->newdir();
	`flash -m $min_overlap -M $max_overlap -I -d $tempdir $input_file`;
	`cp $tempdir/out.notCombined.fastq $paired_output`;
}

