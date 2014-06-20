#!/usr/bin/perl
# April 2014 Takes output of RE blast the output

use strict;
use File::Temp ();
use Bio::SearchIO;
use Getopt::Long;

my $PROPDIV = 10; # when reporting global ouput number of divisions in satellite proportion

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
my @sat_reads; # number of reads associated with satellites in increments of $PROPDIV division of for mininmal proportion of contigs repeats to call a cluster a satellite

### go through clusters

# setup output files
my $percluster_repeat_output = $outputdirectory . "\/" . "satellite_proportion_per_cluster.txt";
open (OUTPUT, ">$percluster_repeat_output") or die "cannot create file $percluster_repeat_output\n";
my $trf_match_assemblies = $outputdirectory . "\/" . "assemblies_with_repeats.fas";
open (OUTPUT2, ">$trf_match_assemblies") or die "cannot create file $trf_match_assemblies\n";

# loop through clusters
foreach my $cluster (@clusterdir) {

	my $contig_file = `find $cluster -name "*.minRD5_sort-GR"`;
	next unless (filesize($contig_file)); #skip the analyis unless the alignment file is not empty

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

	# run trf
	chomp $contig_file;
	`trf $contig_file 1 1 2 80 5 200 2000 -d -h`; # parameters from paper Melters, Daniël P, Keith R Bradnam, Hugh A Young, Natalie Telis, Michael R May, J Ruby, Robert Sebra, et al. “Comparative Analysis of Tandem Repeats from Hundreds of Species Reveals Unique Insights into Centromere Evolution.” Genome Biology 14, no. 1 (2013): R10. doi:10.1186/gb-2013-14-1-r10.
	my @data = split("/", $contig_file);
	my $ouputfile = $data[-1] . ".1.1.2.80.5.200.2000.dat";
	
	my $total_contig_length;
	my $total_repeated_length;
	
	# process trf output
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
			my $print_assembly = 0; # boolean, turns true if this a contig has at least one trf hit then it contains the line to print
			while ($line2 =~ /^\d+\s\d+\s/) {
				my @data = split(" ", $line2);
				my $len = $data[1] - $data[0];
				if ($len > $longest_repeated_length) {
					$longest_repeated_length = $len;
				}
				$numrepeats++;
				$line2 = <INPUT>;

				# updating the line to output
				my $assembly_sequence = $data[-1];
				$print_assembly = ">$title\n" . "$assembly_sequence\n";
			}
			$total_repeated_length += $longest_repeated_length;
			<INPUT>;
			if ($print_assembly) {
				print OUTPUT2 "$print_assembly";
			}
		}
		else {
			die "cannot read line $line";
		}
	}
	
	my $repeatratio = $total_repeated_length/$total_contig_length;
	print OUTPUT "$cluster_shortname\t$repeatratio\n";
	`mv $ouputfile $clusteroutput`;

	for (my $i=0; $i<=$PROPDIV; $i++) {
		if ($repeatratio >= ($i/$PROPDIV)) {
			$sat_reads[$i] += $totalreads_cluster;
		}
	}	
}
close OUTPUT;
close OUTPUT2;

# give global output
open (OUTPUT, ">$outputdirectory/summary_satellites.txt") or die;
print OUTPUT "Column 1: proportion of alignments that have TRF find repeats in to call the whole cluster a satellite\n";
print OUTPUT "Column 2: proportion of total reads in this analysis ($totalreads) that are in clusters identified as satellites\n\n";
for (my $i=0; $i <= scalar @sat_reads; $i++) {
	my $sat_prop = $sat_reads[$i]/$totalreads;
	print OUTPUT"$i\t$sat_prop\n";
}
close OUTPUT;
#print STDERR "proportion of satellite reads $sat_prop\n";


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

# get the size of a file in bytes
sub filesize {
	(my $filename) = @_;
	my $cmdoutput = `ls -l $filename`;
	my @data = split(" ", $cmdoutput);
	return($data[4]);
}
