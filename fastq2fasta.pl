#!/usr/bin/perl
use strict;

use strict;
require File::Temp;
use File::Temp ();
use File::Basename;
use Getopt::Long;
use File::Path;
use List::Util 'shuffle';

my $inputfile; 
my $outputfile;
my $seqlen;
my $numberlines;
my $numprinted=0; # number of sequences printed
my $randomize=0; # pick random lines in the fastq file

#####read and check the inputs
GetOptions(
	'in:s'   => \$inputfile,
	'o:s'	=> \$outputfile,
	'l:s'	=> \$seqlen,
	'n:s'   => \$numberlines,
	'r:s'	=> \$randomize
);

unless ($inputfile and $outputfile) {
	die ("usage: perl fastq2fasta -in <REQUIRED fastqfile, paired> -o <REQUIRED ouptput file, default stdout> -l <exact length of output seauences> -n <number of sequences to output from start of file> -r <number of lines to output picked at random from fastq");
}
if ($randomize and $numberlines) {
	die ("cannot select both -n and -r");
}
unless ($numberlines) {
	$numberlines= 9**9**9; # assigning infinity
}

my %seqstoprint; #holds number of line to print as key
if ($randomize) {
	my $wc = `wc -l $inputfile`;	
	my @data = split(" ", $wc);
	my $fastqseqs = $data[0]/8; # number of sequences in pairs the fastq file
	my @seqnum; #array holding sequence numbers
	for (my $i=0; $i<=$fastqseqs; $i++) {
		push @seqnum, $i;
	}
	my @randseqnum = shuffle @seqnum; # randomize array
	for (my $i=0; $i<$randomize; $i++) {
		$seqstoprint{$randseqnum[$i]}=0;
	}
}

open OUTPUT, ">$outputfile" or die "cannot open $outputfile for writting\n";

###### produce data
my $seqnum=-1;
open (INPUT, $inputfile) or die "cannot open input file $inputfile\n";
while ((my $l1 = <INPUT>) and ($numprinted <= $numberlines)) {
	$seqnum++;
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
	my $skip=0; # if randomize this be used to skip lines not in the array
	if ($randomize) {
		unless (exists $seqstoprint{$seqnum}) {
			$skip=1;
		}
	}
	next if $skip;
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
		$numprinted += 4;
	}	
}
close OUTPUT;
