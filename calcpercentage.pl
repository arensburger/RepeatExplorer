#!/usr/bin/perl
# calculates the percentage TE occupancy
use strict;
use Text::CSV;
use Getopt::Long;

my $manual_annotation_file; # manual annotation file csv
my $REdirectory; # Repeat Explorer output directory name
my $outputfile; # ouput file
GetOptions(
	'm:s'   => \$manual_annotation_file,
	'r:s'	=> \$REdirectory,
	'o:s'	=> \$outputfile,
);
unless ($manual_annotation_file and $REdirectory and $outputfile) {
	print "usage: perl calcpercentage\n\t-m <REQUIRED, manual annotation file as a CSV file>\n\t-r <REQUIRED, Repeat Explorer output directory>\n\t-o <REQUIRED, output file name>\n";
}

### Global variables ###
my %readnumbers; # short cluster name as key number of reads per cluster as value
my %categorienumbers; # repeat category names as key and cluster read abundance as value

## read the number of reads in each cluster
#my @clusterdir = split (" ", `find $REdirectory/seqClust/clustering/clusters/ -name "dir_*"`);
#foreach my $cluster (@clusterdir) {
#	my $totalreads_cluster = `grep -v ">" $cluster/reads.fas -c`;
#	print "$cluster\t$totalreads_cluster\n";
#}
#exit;

## read the manual annotation and populate the %categorienumbers hash
my $csv = Text::CSV->new({ sep_char => ',' }); # create the CSV object
open(my $data, '<', $manual_annotation_file) or die "Could not open '$manual_annotation_file' $!\n";
while (my $line = <$data>) {
	chomp $line;
	if ($csv->parse($line)) {
		my @fields = $csv->fields();
		my $category = $fields[0];
		if (length $fields[1] > 0) {
			$category .= "/" . $fields[1];
		}
		$categorienumbers{$category} += 1;
  	} 
	else {
      		warn "Line could not be parsed: $line\n";
  	}
}

exit;

