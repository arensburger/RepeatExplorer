#!/usr/bin/perl

use strict;

my $i = 0;
open (INPUT, $ARGV[0]) or die;
while (my $line = <INPUT>) {
	$i++;
	unless ($line =~ /^>\S+/) {
		print "line $i is not a title as expected\n";
		exit;
	}
	$i++;
	$line = <INPUT>;
	if ($line =~ /^\s+$/) {
		print "line $i is blank\n";
		exit;
	}
	if ($line =~ /^>/) {
		print "line $i is preceeded by another fasta header\n";
		exit;
	}

}
	

