#!/usr/bin/perl
use strict;

use strict;
require File::Temp;
use File::Temp ();
use File::Basename;
use Getopt::Long;
use File::Path;

my $inputfile; 
my $outputfile;
my $seqlen;
my $numberlines;
my $numprinted=0; # number of sequences printed

#####read and check the inputs
GetOptions(
	'in:s'   => \$inputfile,
	'o:s'	=> \$outputfile,
	'l:s'	=> \$seqlen,
	'n:s'   => \$numberlines,
);

unless ($inputfile and $outputfile) {
	die ("usage: perl fastq2fasta -in <REQUIRED fastqfile, paired> -o <REQUIRED ouptput file, default stdout> -l <exact length of output seauences> -n <number of sequences to output>");
}
unless ($numberlines) {
	$numberlines= 9**9**9; # assigning infinity
}
open OUTPUT, ">$outputfile" or die "cannot open $outputfile for writting\n";

###### produce data
open (INPUT, $inputfile) or die "cannot open input file $inputfile\n";
while ((my $l1 = <INPUT>) and ($numprinted <= $numberlines)) {
	$l1 =~ s/(@)|(\s)//g; 
	my $l2 = <INPUT>;
	 <INPUT>;
	 <INPUT>;
	my $l3 =  <INPUT>;
	$l3 =~ s/(@)|(\s)//g; 
	my $l4 = <INPUT>;
	<INPUT>;
	<INPUT>;
	chomp $l2;
	chomp $l4;
	
	my $printok=0; # boolean for printing this or not
	if ($l1 eq $l3) { 
		if ($seqlen) {
			if (((length $l2) == $seqlen) and ((length $l4) == $seqlen)) {
				$printok=1;	
			}
		}
		else {
			$printok=1;
		}
	}
	else {
		die "sequences are not the same\n$l1\n$l3\n";
	}

	if ($printok) {
		print OUTPUT ">$l1", '\1', "\n";
		print OUTPUT "$l2\n";
		print OUTPUT ">$l3", '\2', "\n";
		print OUTPUT "$l4\n";	
		$numprinted++;
	}	
}
close OUTPUT;
