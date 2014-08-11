#!/usr/bin/perl
use strict;
use Getopt::Long;

my %config;
GetOptions(\%config,
	'redirectory=s',
	'out=s',
);

##> Print USAGE if --help
if ($config{help}) {printUsage(1);}

##> Check if no mandatory parameter is missing and set defaults
if (!exists $config{redirectory})         {printUsage();}
if (!exists $config{out}) {$config{out} = "out";}

##> Open output file
open(OUTPUT, ">$config{out}") or die "cannot open output file $config{out}\n";

# read the number of reads in each cluster
my @clusterdir = split (" ", `find $config{redirectory}/seqClust/clustering/clusters/ -name "dir_*"`);
foreach my $cluster (@clusterdir) {
	my $totalreads_cluster = `grep -v ">" $cluster/reads.fas -c`;
	my $shortclustername = shortname($cluster);
	print OUTPUT "$shortclustername\t$totalreads_cluster";
}
close OUTPUT;

########### subroutines #################################################
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

###########################################################################
sub printUsage{

print STDOUT "DESCRIPTION: This program takes a RepeatExplorer output and reports the number of reads per cluster\n";
print STDOUT "USAGE : RE_readspercluster.pl -f \"reads in fasta format\" -b \"busco file\" -o \"output file\"
Options : 
    -r | redirectory	RepeatExplorer directory (Mandatory)
    -o | out   		Name of ouput file (default \"out\")
    -h | help		You already know it ...\n\n";
    exit;
}
exit;

