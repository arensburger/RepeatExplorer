#!/usr/bin/perl
# June 2014 Takes RE output directory and compares all the assembly files to pFam

use strict;
use File::Temp ();
use Bio::SearchIO;
use Getopt::Long;

my $REP_PROP = 0.5; # minimal proportion of contigs with repeats to call cluster a satellite
my $CPU=8;

### read and check the inputs
my $datadir; # directory with RE output
my $outputdirectory; 

GetOptions(
	'in:s'   => \$datadir,
	'o:s'	=> \$outputdirectory,
);
unless ($datadir and $outputdirectory) {
	die ("usage: perl RE-runPfam -in <RE output directory, REQUIRED> -o <output directory, REQUIRED>\n");
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


### go through clusters
foreach my $cluster (@clusterdir) {
	#get the cluster short name
	my $cluster_shortname = shortname($cluster);
	
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

	#run pfam
	my $translate_output = File::Temp->new( UNLINK => 1, SUFFIX => '.fa' );
	my $pfamoutput =  $clusteroutput . "\/" . "pfam_output.txt";

	`transeq $contig_file $translate_output -frame=6 2>&1`; # translate the DNA into 6 frames
	`~/bin/PfamScan/pfam_scan.pl -fasta $translate_output -dir ~/db/Pfam -pfamB -cpu 8 -outfile $pfamoutput`;

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
	print "$cluster_shortname\t$numhits\t$hits\n"; 
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
