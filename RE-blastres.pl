#!/usr/bin/perl
# April 2014 Takes output of RE blast the output

use strict;
use File::Temp ();
use Bio::SearchIO;
use Getopt::Long;

my $READRM_CUTOFF=0.25; # minimum prortion of RM reads of one class necessary to call the category of a cluster
my $NUMBEROFCPUS = 8;
my $REPBASEDATABASE = "~/db/RepBase/all"; # formated for use with blastn
my $NTDATABASE = "~/db/nt/nt"; # formated for use with blastn

### read and check the inputs
my $datadir; # directory with RE output
my $outputdirectory; 

GetOptions(
	'i:s'   => \$datadir,
	'o:s'	=> \$outputdirectory
);
unless ($datadir and $outputdirectory) {
	die ("usage: perl RE-blastres -i <RE output directory, REQUIRED> -o <output directory, REQUIRED>");
}
`mkdir $outputdirectory`;
if ($?){
	die;
} 

### setup global variables
my @clusterdir = split (" ", `find $datadir/seqClust/clustering/clusters/ -name "dir_*"`);
my %clusterevidence; # cluster name as key, [0] category, [1] number of reads in cluster [2] evidendence of classification based reads [3] evidence based on contigs [4] number of hits based on tblastx of RM


### Step 1 determine which clusters can be assigned to a category based on reads RM 
foreach my $cluster (@clusterdir) {
	#get the cluster short name
	my $cluster_shortname = shortname($cluster);
	
	#create cluster data output directory
	my $clusteroutput =  $outputdirectory . "\/" . $cluster_shortname;
	`mkdir $clusteroutput`;
	if ($?){
		die;
	} 

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
	close INPUT;

	foreach my $read (keys %readclass) {
		my @data2 = split (" ", $readclass{$read});
		foreach my $classcategory (@data2) {
			$class{$classcategory} += 1;
		}
	}

	# report the proportions
	my $first=1; # value goes to 0 after the first line is read
	open (OUTPUT, ">$clusteroutput/RM-read-proportions.txt") or die "cannot create file $clusteroutput/RM-read-proportions.txt";
	foreach my $category (sort {$class{$b} <=> $class{$a}} (keys %class)) {
		my $proportion = $class{$category}/$totalreads_cluster;
		if (($first) and ($proportion >= $READRM_CUTOFF)){
			$clusterevidence{$cluster}[0] = $category;
			$clusterevidence{$cluster}[1] = $totalreads_cluster;
			$clusterevidence{$cluster}[2] = $proportion;
		}
		print OUTPUT "$category\t$proportion\n";
		$first = 0;
	}
	close OUTPUT;
	unless (exists $clusterevidence{$cluster}[0]) { # go here if a classification could not be established
		$clusterevidence{$cluster}[0] = 0;
		$clusterevidence{$cluster}[1] = $totalreads_cluster;
	}
}

### Step 2 determine which clusters can be assigned to a category based on contig similarity to known TEs
foreach my $cluster (keys %clusterevidence) { 
	my $clusteroutput =  $outputdirectory . "\/" . shortname($cluster);
	if ($clusterevidence{$cluster}[0] eq 0) {
		my $numcontigs; #total number of contigs for this cluster

		### add all the contigs into a temporary file ###
		my $rm_blastinput = File::Temp->new( UNLINK => 1, SUFFIX => '.fa' );
		open (OUTPUT, ">$rm_blastinput") or die "cannot open output file $rm_blastinput";
		open (INPUT, `find $cluster -name *.minRD5_sort-GR`) or die "cannot open file `find $cluster -name *.minRD5_sort-length`";
		
		## read the data from the file
		# process first line	
		my $line = <INPUT>;
		chomp $line;
		print OUTPUT "$line\n";
	
		# process remaining lines
		my $sequence;
		while (my $line = <INPUT>) {
			if ($line =~ />/) {
				print OUTPUT "$sequence\n";
				$numcontigs++;
				$sequence = "";
				chomp $line;
				print OUTPUT "$line\n";
			}
			else {
				chomp $line;
				$sequence .= $line;
			}
		}
		print OUTPUT "$sequence\n";
		$numcontigs++;
		close INPUT;
		close OUTPUT;

		# perform the blast on RM using nucleotide blast
		my $rm_blastoutput = "$clusteroutput/contig-rm-blast.txt";
		`blastn -query $rm_blastinput -db $REPBASEDATABASE -out $rm_blastoutput -outfmt 6 -evalue 0.1 -num_alignments 1 -num_threads $NUMBEROFCPUS 2>/dev/null`;
		# parse blast output
		my %contigs; # list of contigs that had match
		my %repeats; # list of repeats that matched the contigs

		my $rmblastres; #results from blast of contigs to repeatmasker
		open (INPUT, $rm_blastoutput) or die "cannot open $rm_blastoutput";
		while (my $line = <INPUT>) {
			my @data = split (" ", $line);
			$contigs{$data[0]} = 0;
			$repeats{$data[1]} = 0;
		}
		my $numhitcontigs = keys %contigs;
		my $numrepeats = keys %repeats;

		if ($numhitcontigs > 0) {
			$rmblastres = "of $numcontigs contigs $numhitcontigs matched repeats: ";
			foreach my $key (keys %repeats) {
				$rmblastres .= "$key ";
			}
			$clusterevidence{$cluster}[0] = "contig RM hit";
			$clusterevidence{$cluster}[3] = $rmblastres;
		}

		# perform the blast on RM using tblastx
		my $tblastxoutput = "$clusteroutput/rmtblastx.xml";
		`tblastx -query $rm_blastinput -db $REPBASEDATABASE -out $tblastxoutput -outfmt 5 -evalue 0.1 -num_alignments 1 -num_threads $NUMBEROFCPUS 2>/dev/null`;
		my $tblastxsimple = "$clusteroutput/rmtblastx.xls";
		open (OUTPUT, ">$tblastxsimple") or die "cannot create file $tblastxsimple\n";
		my $counter = 0;
		my $searchin = new Bio::SearchIO( -tempfile => 1,
				 		 -format => 'blastxml',
				 		 -file   => $tblastxoutput);
		while( my $result = $searchin->next_result ) {
			my $query_name = $result->query_description;
			while (my $hit = $result->next_hit) {
				my $hit_name = $hit->name;
				my $hsp = $hit->next_hsp;
				my $evalue = $hsp->evalue;
				my $description = $hit->description;
				my $start_query = $hsp->start('query');
				my $end_query = $hsp->end('query');
				$counter++;
				print OUTPUT "$query_name\t$hit_name\t$evalue\t$start_query\t$end_query\n";
			}

		}
		if ($counter) {
			$clusterevidence{$cluster}[4] = $counter;
		}
		close OUTPUT;

	}
}

### Step 2-2 run blast on remaining clusters
foreach my $cluster (keys %clusterevidence) { 
	my $clusteroutput =  $outputdirectory . "\/" . shortname($cluster);
	if (($clusterevidence{$cluster}[0] eq 0) or ($clusterevidence{$cluster}[0] eq "contig RM hit")) {
		my $blastoutput = "$clusteroutput/contig-nt-blast.xml";
		my $blastinput = `find $cluster -name *.minRD5_sort-GR`;
		chomp $blastinput;
		`blastn -query $blastinput -db $NTDATABASE -out $blastoutput -evalue 0.001 -num_threads $NUMBEROFCPUS -outfmt 5 2>/dev/null`;	

		# make short output version
		if (-e $blastoutput) { #test if the blast created a file or not
			my $counter = 0;
			my $summaryoutput = "$clusteroutput/contig-nt-blast.xls";
			open (OUTPUT, ">$summaryoutput") or die "cannot create file $summaryoutput\n";
			print OUTPUT "name\tevalue\thit name\tquery start\tquery end\n";
			my $searchin = new Bio::SearchIO( -tempfile => 1,
					 		 -format => 'blastxml',
					 		 -file   => $blastoutput);
			while( my $result = $searchin->next_result ) {
				my $query_name = $result->query_description;
				while (my $hit = $result->next_hit) {
					my $hit_name = $hit->name;
					my $hsp = $hit->next_hsp;
					my $evalue = $hsp->evalue;
					my $description = $hit->description;
					my $start_query = $hsp->start('query');
					my $end_query = $hsp->end('query');
					$counter++;
					print OUTPUT "$query_name\t$evalue\t$description\t$start_query\t$end_query\n";
				}
	
			}
			if ($counter) {
				$clusterevidence{$cluster}[5] = $counter;
			}
			close OUTPUT;
		}
	}
}


### Step 3 report data into excel database
my $summaryfile = $outputdirectory . "/" . "summary.xls";
open (OUTPUT, ">$summaryfile") or die "cannot create file $summaryfile\n";
foreach my $cluster (sort {$clusterevidence{$b}[1] <=> $clusterevidence{$a}[1]} (keys %clusterevidence)) { #sort by number of reads
	if ($clusterevidence{$cluster}[0] eq "contig RM hit") {
		my $shortname = shortname($cluster);
		print OUTPUT "$shortname\t$clusterevidence{$cluster}[1]\t\t$clusterevidence{$cluster}[4]\t$clusterevidence{$cluster}[5]\t$clusterevidence{$cluster}[3]\n";
	}
	elsif ($clusterevidence{$cluster}[0] eq 0) {
		my $shortname = shortname($cluster);
		print OUTPUT "$shortname\t$clusterevidence{$cluster}[1]\t\t$clusterevidence{$cluster}[4]\t$clusterevidence{$cluster}[5]\n";
	}
	else {
		my $shortname = shortname($cluster);
		print OUTPUT "$shortname\t$clusterevidence{$cluster}[1]\t$clusterevidence{$cluster}[0]\n";
	}
}
close OUTPUT;


### step 4 report overall stats
my $totreads = (split ' ', `grep "Formatted" $datadir/seqClust/sequences/formatdb.log`)[1]; # total reads in analyis
my $cluster_abund; # number of reads in clusters
my $cluster_reads; # number of reads associated with clusters that can be categorised based on read RM
my $cluster_contigs; # number of reads associated with clusters that can be categorised based on contig blast
my $num_clusters; # number of clusters
my $num_clusters_reads; # number of clusters that can be categorized by reads RM
my $num_clusters_contig; # number of clusters that can be categorized by contigs blast

foreach my $cluster (keys %clusterevidence) {
	$cluster_abund += $clusterevidence{$cluster}[1];
	$num_clusters++;
	if ($clusterevidence{$cluster}[0] eq "contig RM hit") {
		$cluster_contigs += $clusterevidence{$cluster}[1];
		$num_clusters_contig++;
	}
	elsif ($clusterevidence{$cluster}[0] eq 0) {
		
	}
	else {
		$cluster_reads += $clusterevidence{$cluster}[1];
		$num_clusters_reads++;
	}
}

my $prop_reads = $cluster_reads/$cluster_abund;
my $prop_contigs = $cluster_contigs/$cluster_abund;
print "Number of reads $totreads, reads in clusters $cluster_abund\n";
print "$prop_reads of cluster reads can be categorized by read RM\n";
print "$prop_contigs of clusters reads might be categorized by contig blast to RM database\n";
print "Summary file is in $summaryfile\n";


#### SUBROUTINES ########

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

# get the short name of the cluster
sub shortname {
	(my $cluster) = @_;
	my $shortname;
	if ($cluster =~ /dir_(CL\d+)/) {
		$shortname = $1;
	}
	else {
		die "cannot parse cluster name: $cluster";
	}
	return ($shortname);
}
