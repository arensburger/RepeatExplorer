#!/usr/bin/perl

use strict;

my $i = 0;
my $ref1;
my $ref2;
open (INPUT, $ARGV[0]) or die;
my $ref1 = <INPUT>;
my $ref2  = <INPUT>;
while (my $l1 = <INPUT>) {
	$i++;
	my $l2 = <INPUT>;
	$i++;

	my $sref1 = substr($ref1, 0, -3);
	my $sl1 = substr($l1, 0, -3);
	if ($sref1 eq $sl1) {
		print "$sref1", '\1', "\n";
		print $ref2;
		print "$sl1", '\2', "\n";
		print $l2
	}
	else {
		$ref1 = $l1;
		$ref2 = $l2;
	}	
}

