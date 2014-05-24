#!/usr/bin/perl

use strict;

my %ele;
open (INPUT, $ARGV[0]) or die;
while (my $line = <INPUT>) {
	my @data = split (" ", $line);
	$ele{$data[2]} += $data[1];
}
my $percent;
foreach my $key (keys %ele) {
	if ($ARGV[1]) {
		$percent = $ele{$key}/$ARGV[1];
	}
	else {
		$percent = "N/A";
	}
	print "$key\t$ele{$key}\t$percent\n";
}
