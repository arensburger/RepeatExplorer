#!/usr/bin/perl
# June 2014 Takes RE output directory and compares all the assembly files to pFam

use strict;
use File::Temp ();
use Bio::SearchIO;
use Getopt::Long;

my $CPU=2;
my $pfamB=0; # set to something else than 0 to run on both pfamA and B databases

### set default database locations ###
my $pfamDB="~/db/Pfam";
my $pfamREPET="~/db/ProfilesBankForREPET/ProfilesBankForREPET_Pfam26.0_GypsyDB.hmm";

### read and check the inputs
my $datadir; # directory with RE output
my $outputdirectory; 
my $file_clusternames;

GetOptions(
	'in:s'   => \$datadir,
	'o:s'	=> \$outputdirectory,
	'b:s'	=> \$pfamB,
	'f:s'	=> \$file_clusternames,
	'pfamDB:s'	=> \$pfamDB,
	'pfamREPET:s'	=> \$pfamREPET,
);
unless ($datadir and $outputdirectory) {
	die ("usage: perl RE-runPfam\n\t-in <RE output directory, REQUIRED>\n\t-o <output directory, REQUIRED>\n\t-b <run against pfamB in addition to pfamA>\n\t-f <file with names of clusters to anlalyse>\n\t-pfamDB <location of pfam database, currently set to $pfamDB>\n\t-pfamREPET <location of REPET currated hmm file, currently set to $pfamREPET>\n");
}
# if list file with list of clusters is provided then populate the names into array
my %clustertocheck; #names of clusters to check (if any)
if ($file_clusternames) {
	open (INPUT, $file_clusternames) or die "cannot open file $file_clusternames\n";
	while (my $line = <INPUT>) {
		chomp $line;
		$clustertocheck{$line} = 1;
	}
	close INPUT;
}

`mkdir -p $outputdirectory`;
if ($?){
	die;
} 

### setup global variables
my @clusterdir = split (" ", `find $datadir/seqClust/clustering/clusters/ -name "dir_*"`);
# get total number of reads
my $line = `grep Formatted $datadir/seqClust/sequences/formatdb.log`;
my @data = split (" ", $line);
my $totalreads = $data[1];
my $sat_reads; # number of reads associated with satellites 
my $summaryfilename = $outputdirectory . "\/" . "pfamSummary.xls";
open (OUTPUT, ">$summaryfilename") or die "cannot make file $summaryfilename\n";

### go through clusters
foreach my $cluster (@clusterdir) {
	#get the cluster short name
	my $cluster_shortname = shortname($cluster);
	
	my $ignorecluster=0; # boolean set to 1 if ignore cluster
	if ($file_clusternames) {
		unless ($clustertocheck{$cluster_shortname}) {
			$ignorecluster=1;
		}
	}
	next if($ignorecluster);

	#number of reads in the cluster
	my $totalreads_cluster = `grep -v ">" $cluster/reads.fas -c`;
	chomp $totalreads_cluster;

	#create cluster data output directory
	my $clusteroutput =  $outputdirectory . "\/" . $cluster_shortname;
	`mkdir $clusteroutput`;
	if ($?){
		die;
	} 
	my $contig_file = `find $cluster -name "*.minRD5_sort-GR"`;
	chomp $contig_file;

	### run the pfams
	my $translate_output = File::Temp->new( UNLINK => 1, SUFFIX => '.fa' );
	my $pfamoutput =  $clusteroutput . "\/" . "pfam_output.txt";
	my $pfamREPEToutput = $clusteroutput . "\/" . "pfamREPET_output.txt";

	`transeq $contig_file $translate_output -frame=6 2>&1`; # translate the DNA into 6 frames

	# standard pfam
	if ($pfamB) {
		`pfam_scan.pl -fasta $translate_output -dir $pfamDB -pfamB -cpu $CPU -outfile $pfamoutput`;
	}
	else {
		`pfam_scan.pl -fasta $translate_output -dir $pfamDB -cpu $CPU -outfile $pfamoutput`;
	}

	# parse the pfam output	
	my $numhits; #number of hits
	my %hmmname; # name of hit at key and frequency as value
	open (INPUT, $pfamoutput) or die "cannot open file $pfamoutput\n";
	while (my $line = <INPUT>) {
		unless (($line =~ /^#/) or (length $line < 2)) {
			my @data = split (" ", $line);
			$hmmname{$data[6]} += 1;
			$numhits++;
		}
	}
	close INPUT;
	my $hits;
	foreach my $key (sort { $hmmname{$b} <=> $hmmname{$a} } (keys %hmmname)) {
		my $proportion = $hmmname{$key}/$numhits;
		my $percent = prop2percent($proportion) . "%";
		$hits .= "$key($percent) ";
	}

	# REPET currated
        `hmmscan --tblout $pfamREPEToutput $pfamREPET $translate_output`;

        # parse the pfam output 
        my $numhits; #number of hits
        my %hmmname; # name of hit at key and frequency as value
        open (INPUT, $pfamREPEToutput) or die "cannot open file $pfamoutput\n";
        while (my $line = <INPUT>) {
                unless (($line =~ /^#/) or (length $line < 2)) {
                        my @data = split (" ", $line);
                        $hmmname{$data[0]} += 1;
                        $numhits++;
                }
        }
        close INPUT;
        my $REPEThits;
        foreach my $key (sort { $hmmname{$b} <=> $hmmname{$a} } (keys %hmmname)) {
                my $proportion = $hmmname{$key}/$numhits;
                my $percent = prop2percent($proportion) . "%";
                $REPEThits .= "$key($percent) ";
        }

	print OUTPUT "$cluster_shortname\t$hits\t$REPEThits\n"; 
}

# get the short name of the cluster
sub shortname {
	(my $cluster) = @_;
	my $shortname;
	if ($cluster =~ /dir_(CL\d+)/) {
		$shortname = $1;
	}
	else {
		die "cannot parse cluster name: $cluster";
	}
	return ($shortname);
}

# convert proportion to percent
sub prop2percent {
	(my $num) = @_;
	my $percent = 100 * $num;
	my $rounded = sprintf "%.2f", $percent;
	return $rounded;
}
