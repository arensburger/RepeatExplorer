#!/usr/bin/perl

### Sept 2014. Takes output of RE run and returns a file listing the clusters to merge 

use strict;
#use diagnostics;
use warnings;
use Getopt::Long;
#use agutils;
#use File::Temp ();
use File::Temp qw(tempdir);

##> Define options
my %config;
GetOptions (\%config,
            'in=s',
	    'cutoff=s',
            'out=s',
            'help');

##> Print USAGE if --help
if ($config{help}) {printUsage(1);}

##> Check if no mandatory parameter is missing and set defaults
if (!exists $config{in})         {printUsage();}
if (!exists $config{cutoff})         {$config{cutoff} = 0.1;}
if (!exists $config{out}) {$config{out} = "merge.txt";}
unless (($config{cutoff} == 0.2) or ($config{cutoff} == 0.1) or ($config{cutoff} == 0.05)) {
	printUsage();
}
open (OUTPUT, ">$config{out}") or die "cannot open output file $config{out}\n";
### Main program ###
my $filename = "$config{in}/summary/merge_cutoff" . $config{cutoff} . "/mate_pairs_merging.html";
open (FILE, "$filename") or die "cannot open file $filename\n";
while (my $line = <FILE>) {
	if ($line =~ /^<p\sclass='character'>CL\d+&nbsp;.+<\/p>/) {
		while ($line =~ /(CL\d+)/g) {
			print OUTPUT "$1\t";
		}
		print OUTPUT "\n";
	}
#	if ($line =~ /Individual\sclusters\s-\sconnection\sabsent\sor\sbelow\scutoff/) {
#		print "hqh\n";
#	}
	last if ($line =~ /Individual\sclusters\s-\sconnection\sabsent\sor\sbelow\scutoff/);
}
close OUTPUT;
close FILE;

###########################################################################
################################ Functions ################################
###########################################################################
sub printError{
    my $string = shift;
    my $exit = shift;
    
    print STDERR $string;
    exit if $exit;
}

###########################################################################
sub printUsage{

print STDOUT "USAGE : RE-merge-clusters.pl -in \"repeat explorer output directory\" -c \"cutoff\" -o \"output file\"
Options : 
    -in| 		Repeat Explore outut directory (Mandatory)
    -c | cutoff	 	Cluster assembly cutoff 0.2, 0.1, or 0.05 (default 0.1)
    -o | out   		Name of ouput file (default \"merge.txt\")
    -h | help		You already know it ...\n\n";
    exit;
}
