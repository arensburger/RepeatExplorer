#!/usr/bin/perl
# April 2014 Takes output of RE blast the output

use strict;
use File::Temp ();
use Bio::SearchIO;
use Getopt::Long;

my $REP_PROP = 0.5; # minimal proportion of contigs with repeats to call cluster a satellite

### read and check the inputs
my $datadir; # directory with RE output
my $outputdirectory; 

GetOptions(
	'in:s'   => \$datadir,
	'o:s'	=> \$outputdirectory,
);
unless ($datadir and $outputdirectory) {
	die ("usage: perl RE-blastres -in <RE output directory, REQUIRED> -o <output directory, REQUIRED>\n");
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
	`trf $contig_file 1 1 2 80 5 200 2000 -d -h`; # parameters from paper Melters, Daniël P, Keith R Bradnam, Hugh A Young, Natalie Telis, Michael R May, J Ruby, Robert Sebra, et al. “Comparative Analysis of Tandem Repeats from Hundreds of Species Reveals Unique Insights into Centromere Evolution.” Genome Biology 14, no. 1 (2013): R10. doi:10.1186/gb-2013-14-1-r10.
	my @data = split("/", $contig_file);
	my $ouputfile = $data[-1] . ".1.1.2.80.5.200.2000.dat";
	
	my $total_contig_length;
	my $total_repeated_length;
	
	open (INPUT, $ouputfile) or die "cannot open file $ouputfile\n";
	<INPUT>;
	<INPUT>;
	<INPUT>;
	<INPUT>;
	<INPUT>;
	<INPUT>;
	<INPUT>;
	<INPUT>;
	while (my $line = <INPUT>) {
		my $title; # name of current sequence
		my $length_contig; # length of current contig
		my $length_rep; # length of current contig covered by repeats
		my $numrepeats = 0; # nubmer of different repeats
		if ($line =~ /Sequence:\s(\S+)\s\((\d+)-/) {
			$title = $1;
			$total_contig_length += $2;
			<INPUT>;
			<INPUT>;
			<INPUT>;
			<INPUT>;
			<INPUT>;
			<INPUT>;
			my $line2 = <INPUT>;
			my $longest_repeated_length;
			while ($line2 =~ /^(\d+)\s(\d+)\s/) {
				my $len = $2 - $1 + 1;
				if ($len > $longest_repeated_length) {
					$longest_repeated_length = $len;
				}
				$numrepeats++;
				$line2 = <INPUT>;
			}
			$total_repeated_length += $longest_repeated_length;
			<INPUT>;
		}
		else {
			die "cannot read line $line";
		}
	}
	
	my $repeatratio = $total_repeated_length/$total_contig_length;
	print "$cluster_shortname\t$repeatratio\n";
	`mv $ouputfile $clusteroutput`;

	if ($repeatratio >= $REP_PROP) {
		$sat_reads += $totalreads_cluster;
	}
}

# give global output
my $sat_prop = $sat_reads/$totalreads;
print STDERR "proportion of satellite reads $sat_prop\n";



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
