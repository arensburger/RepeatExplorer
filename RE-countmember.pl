#!/usr/bin/perl
# April 2014 Takes output of RE and counts cluster membership

use strict;

my $CLUSTERTHRESHOLD = 0.05; # proportion of reads of one class to call a cluster one element

my $datadir = $ARGV[0]; # directory with RE output
my @clusterdir = split (" ", `find $datadir/seqClust/clustering/clusters/ -name "dir_*"`);
my %membership; # class name as key and number of reads as value
my $totreads = (split ' ', `grep "Formatted" $datadir/seqClust/sequences/formatdb.log`)[1]; # total reads in analyiss

#determine get basic data for all clusters
foreach my $cluster (@clusterdir) {

	#number of reads in the cluster
	my $totalreads_cluster = `grep -v ">" $cluster/reads.fas -c`;
	$totalreads_cluster =~ s/\s//g;

	#parse the RM output
	my %class; # repeat class ID as key and abundance as value

	# count the abundance of classes
	open (INPUT, "$cluster/RM-reads.fas.out") or die "cannot open file $cluster/RM-reads.fas.out\n";
	<INPUT>;
	<INPUT>;
	<INPUT>;
	while (my $line = <INPUT>) {
		my @data = split (" ", $line);
		$class{classformat($data[10])} += 1;
	}

	# get the large value in %class, from http://www.perlmonks.org/bare/?node_id=539366
	my @key_list = keys %class;
	my $max_key = pop @key_list;
	foreach (@key_list) {
    		$max_key = $_ if $class{$_} > $class{$max_key};
	}

	# determine if the reads in this cluster should be added to the total
	if (($class{$max_key}/$totalreads_cluster) >= $CLUSTERTHRESHOLD) {
		$membership{$max_key} += $totalreads_cluster;
	}
}

foreach my $key (keys %membership) {
	my $percentage = ($membership{$key}/$totreads)*100;
	print "$key\t$membership{$key}\t$percentage", '%', "\n";
}

#format the class name
sub classformat {
	my ($text) = @_;
	my $formated_text;

	if ($text =~ /(\S+)\//) {
		$formated_text = $1;		
	}
	else {
		$formated_text = $text;
	}
	
	return $formated_text;
}
