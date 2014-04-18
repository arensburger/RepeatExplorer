#!/usr/bin/perl
# April 2014 Takes output of RE blast the output

use strict;
use File::Temp ();
use Bio::SearchIO;

my $TOPREADS=3; # number of best reads to report
my $CLUSTERTHRESHOLD = 0.05; # proportion of reads of one class to call a cluster one element

my $datadir = $ARGV[0]; # directory with RE output
my @clusterdir = split (" ", `find $datadir/seqClust/clustering/clusters/ -name "dir_*"`);
my %membership; # class name as key and number of reads as value
my $totreads = (split ' ', `grep "Formatted" $datadir/seqClust/sequences/formatdb.log`)[1]; # total reads in analyiss

#determine get basic data for all clusters
foreach my $cluster (@clusterdir) {

	# get the cluster name
	my $clustername;
	if ($cluster =~ /dir_(CL\d+)/) {
		$clustername = $1;
	}
	else {
		die "cannot get cluster name from\n$cluster\n";
	}

	### count the reads from the RM output ###

	# data collected
	my @top_rm_read_classes; # name of the top RM read classes
	my @top_rm_read_numbers; # number of RM reads for each top class

	#number of reads in the cluster
	my $totalreads_cluster = `grep -v ">" $cluster/reads.fas -c`;
	$totalreads_cluster =~ s/\s//g;

	# count the abundance of classes, need to classify each read individually because the same read can have multiple classes
	my %readclass; #associates each read to one or more classes
	my %class; # repeat class ID as key and abundance as value, these are totals
	open (INPUT, "$cluster/RM-reads.fas.out") or die "cannot open file $cluster/RM-reads.fas.out\n";
	<INPUT>;
	<INPUT>;
	<INPUT>;
	while (my $line = <INPUT>) {
		my @data = split (" ", $line);
		my $currentclass = classformat($data[10]);
		unless ($readclass{$data[4]} =~ /$currentclass/) {
			$readclass{$data[4]} .= "$currentclass ";
		}
		$class{$currentclass} = 0; # record that this class exists
	}
	foreach my $read (keys %readclass) {
		my @data2 = split (" ", $readclass{$read});
		foreach my $classcategory (@data2) {
			$class{$classcategory} += 1;
		}
	}
	my $i=1;
	foreach my $key (sort {$class{$b} <=> $class{$a}} (keys %class)) {
		if ($i<=$TOPREADS) {
			push @top_rm_read_classes, $key;
			my $percentage = sprintf "%.2f", 100*$class{$key}/$totalreads_cluster;
			push @top_rm_read_numbers, $percentage;
		}
		$i++;
	}
	close INPUT;

	### perform blast on best contigs ###
	my @contigname; # holds name of contig in order
	my @contigseq; # holds sequences of contig in order
	open (INPUT, `find $cluster -name *.minRD5`) or die "cannot open file `find $cluster -name *.minRD5`";
	
	## read the data from the file
	# process first line	
	my $line = <INPUT>;
	chomp $line;
	push @contigname, $line;
	
	# process remaining lines
	my $sequence;
	while (my $line = <INPUT>) {
		if ($line =~ />/) {
			push @contigseq, $sequence;
			$sequence = "";
			chomp $line;
			push @contigname, $line;
		}
		else {
			chomp $line;
			$sequence .= $line;
		}
	}
	push @contigseq, $sequence;
	close INPUT;

	## perform the blasts

	# data collected
	my @top_blast_hit; # id and evalue of top bast hit for biggest contigs for RM
	my @top_blastnt_hit; # id and evalue of top bast hit for biggest contigs for nt

	for (my $i=0; $i < $TOPREADS; $i++) {
		my $rm_blastinput = File::Temp->new( UNLINK => 1, SUFFIX => '.fa' );
		my $rm_blastoutput = File::Temp->new( UNLINK => 1, SUFFIX => '.txt' );

		open (OUTPUT, ">$rm_blastinput") or die "cannot open temp file $rm_blastinput";
		print OUTPUT "$contigname[$i]\n$contigseq[$i]\n";
		close OUTPUT;
		
		# RM blast
		`blastn -query $rm_blastinput -db ~/db/RepBase/all -out $rm_blastoutput -outfmt 6 -evalue 0.1 -num_alignments 1 2>/dev/null`;
		my $blastfail = 1 if $?; # capture failure of the bast
		
		unless ($blastfail) {
			my $topline = `head -n 1 $rm_blastoutput`;
			if ($topline eq '') {
				push @top_blast_hit, "no_hit\t";
			}
			else {
				my @data3 = split " ", $topline;
				push @top_blast_hit, "$data3[1]\t$data3[10]"
			}		
		}

		# nt blast
		`blastn -query $rm_blastinput -db ~/db/nt/nt -out $rm_blastoutput -evalue 0.001 -outfmt 5 -num_alignments 1 2>/dev/null`;
 		my $blastfail = 1 if $?; # capture failure of the bast
		
		unless ($blastfail) {
			my $searchin = new Bio::SearchIO( -tempfile => 1,
				  -format => 'blastxml',
				  -file   => $rm_blastoutput);
			while( my $result = $searchin->next_result ) {
				while (my $hit = $result->next_hit) {
					my $hit_name = $hit->name;
					my $hsp = $hit->next_hsp;
					my $evalue = $hsp->evalue;
					my $description = $hit->description;
					push @top_blastnt_hit, "$description $hit_name\t$evalue"
#					print "$clustername\t$hit_name\t$evalue\t$description\n";
				}

			}	
		}
	}

	## output the results
	print "$clustername";

	# print the RM on reads
	for (my $i=0; $i<$TOPREADS; $i++) {
		my $readclass; # class based on reads
		if ($top_rm_read_classes[$i] eq '') {
			$readclass = "NA";
		}
		else {
			$readclass = $top_rm_read_classes[$i];
		}

		my $readpercentage; # percent of total reads for this class
		if ($top_rm_read_numbers[$i] eq '') {
			$readpercentage = "NA";
		}
		else {
			$readpercentage = $top_rm_read_numbers[$i] . '%';
		}

		print "\t$readclass\t$readpercentage";
	}

	# print blast using RM database
	for (my $i=0; $i<$TOPREADS; $i++) {
		if ($top_blast_hit[$i] eq '') {
			print "\t\t";
		}
		else {
			print "\t$top_blast_hit[$i]";
		}
	}

	# print blast using NT database
	for (my $i=0; $i<$TOPREADS; $i++) {
		if ($top_blastnt_hit[$i] eq '') {
			print "\t\t";
		}
		else {
			print "\t$top_blastnt_hit[$i]";
		}
	}

	print "\t$totalreads_cluster\n";
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
