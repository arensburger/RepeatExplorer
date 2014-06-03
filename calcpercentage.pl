#!/usr/bin/perl
# calculates the percentage TE occupancy
use strict;

my %occupancy;
my $seqinclusters; # total number of sequences in clusters

open (INPUT, $ARGV[0]) or die;
while (my $line = <INPUT>) {
	if ($line =~ /\S+\t(\d+)\t(\S*)/) {
		my $abundance = $1;
		my $category = $2;

		if ($category eq "") {
			$occupancy{"not-classified"} += $abundance;
		}
		else {
			$occupancy{$category} += $abundance;
		}
		$seqinclusters += $abundance;
	}
	else {
		warn ("WARNING, not reading line $line");
	}
}

foreach my $cat (keys %occupancy) {
	print "$cat\t$occupancy{$cat}";
	if (exists $ARGV[1]) {
		my $percentage = $occupancy{$cat}/$ARGV[1];
		print "\t$percentage\n";
	}
	else {
		print "\n";
	}
}
if (exists $ARGV[1]) {
	my $percentageincluster = $seqinclusters/$ARGV[1];
	print "Num sequences in clusters\t$seqinclusters\t$percentageincluster\n";
}	
