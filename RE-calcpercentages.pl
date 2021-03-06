#!/usr/bin/perl
# calculates the percentage TE occupancy
use strict;
use Text::CSV;
use Getopt::Long;
use String::ShellQuote;

my $manual_annotation_file; # manual annotation file csv
my $REdirectory; # Repeat Explorer output directory name
my $outputfile; # ouput file
my $clusterabundancefile; # alternative to the REdirectory, this could be the output of RE-readspercluster.pl
my $totreads; # total reads in the analyis

GetOptions(
	'm:s'   => \$manual_annotation_file,
	'r:s'	=> \$REdirectory,
	'f:s'	=> \$clusterabundancefile,
	't:s'	=> \$totreads,
	'o:s'	=> \$outputfile,
);
unless (($manual_annotation_file and $outputfile) and ($REdirectory or $clusterabundancefile))  {
	print "usage: perl calcpercentage\n\t-m <REQUIRED, manual annotation file as a CSV file>\n\t-r <Repeat Explorer output directory>\n\t-f <tab delimited file with cluster abundance, alternative to  RE output directory>\n\t-t <total number of reads>\n\t-o <REQUIRED, output file name>\n\nNote: cannot use spaces in directory names because is screws up the find command, something that needs to be fixed\n";
}

### Global variables ###
my %readnumbers; # short cluster name as key number of reads per cluster as value
my %categorienumbers; # repeat category names as key and cluster read abundance as value
my %clusternumbers; # repeat category names as key and number of clusters as value

# read the number of reads in each cluster
if ($REdirectory) {
	my @clusterdir = split (" ", `find $REdirectory/seqClust/clustering/clusters/ -name "dir_*"`);
	foreach my $cluster (@clusterdir) {
		my $totalreads_cluster = `grep -v ">" $cluster/reads.fas -c`;
		$readnumbers{shortname($cluster)} = $totalreads_cluster;
	}
	$totreads = (split ' ', `grep "Formatted" $REdirectory/seqClust/sequences/formatdb.log`)[1]; # total reads in analyis
}
elsif ($clusterabundancefile) {
	open (INPUT, "$clusterabundancefile") or die "cannot open cluster abundance file $clusterabundancefile\n";
	while (my $line = <INPUT>) {
		my @data = split (" ", $line);
		$readnumbers{$data[0]} = $data[1];
	}
}
else {
	die "ERROR, no cluster abundance provided\n";
}

# read the manual annotation and populate the %categorienumbers hash
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

		if (exists $readnumbers{$fields[2]}) {
			if ($category eq "") {
				$category = "Not classified";
			}
			$categorienumbers{$category} += $readnumbers{$fields[2]};
		}
		else {
			warn "Error, cannot find cluster name $fields[2] from manual annotation in directory structure of $REdirectory; ignoring it\n";
		}
		$clusternumbers{$category} += 1;
  	} 
	else {
      		warn "Line could not be parsed: $line\n";
  	}
}

# output the results
open (OUTPUT, ">$outputfile") or die "cannot open ouput file $outputfile\n";
print OUTPUT "Total number of reads $totreads\n";
print OUTPUT "Category\tProportion of reads\tNumber of cluster\n\n";
foreach my $category (keys %categorienumbers) {
	my $proportion = $categorienumbers{$category}/$totreads;
	print OUTPUT "$category\t$proportion\t$clusternumbers{$category}\n";
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
exit;

