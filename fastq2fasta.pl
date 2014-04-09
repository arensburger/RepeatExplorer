#!/usr/bin/perl
use strict;

open (INPUT, $ARGV[0]) or die;
while (my $l1 = <INPUT>) {
	$l1 =~ s/(@)|(\s)//g; 
	my $l2 = <INPUT>;
	 <INPUT>;
	 <INPUT>;
	my $l3 =  <INPUT>;
	$l3 =~ s/(@)|(\s)//g; 
	my $l4 = <INPUT>;
	<INPUT>;
	<INPUT>;
	print ">$l1", '\1', "\n";
	print $l2;
	print ">$l3", '\2', "\n";
	print $l4;
}
