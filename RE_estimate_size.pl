#!/usr/bin/perl -w

# July 2014 Using BUSCOs as a reference this script takes an set of reads as input and estimates coverage of BUSCOs thus estimating genome percentage

use strict;
use diagnostics;
#use warnings;
use Getopt::Long;
use agutils;
#use File::Temp ();
use File::Temp qw(tempdir);

##> Define options
my %config;
GetOptions (\%config,
            'fasta=s',
	    'busco=s',
	    'threads=s',
            'out=s',
            'verbose',
            'help');

##> Print USAGE if --help
if ($config{help}) {printUsage(1);}

##> Check if no mandatory parameter is missing and set defaults
if (!exists $config{fasta})         {printUsage();}
if (!exists $config{busco})         {$config{busco} = "/home/parensburge/db/BUSCO/BUSCO_Vertebrata.fas";}
if (!exists $config{threads})         {$config{threads} = 12;}
if (!exists $config{out}) {$config{out} = "out";}

## create blast database for BUSCO
my $tempdir = tempdir( CLEANUP => 1 );
print "creating blast database BUSCO using $config{busco} in $tempdir\n" if $config{verbose};
my $cmdoutput = `makeblastdb -in $config{busco} -dbtype prot -out $tempdir/BUSCO`;

## run the blast
my $blast_output = File::Temp->new( UNLINK => 1, SUFFIX => '.blast' ); # temporary file with results blast
print "running blastx\n" if $config{verbose};
$cmdoutput = `blastx -db $tempdir/BUSCO -query $config{fasta} -evalue 1 -num_alignments 1 -num_threads $config{threads} -out $blast_output -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen slen"`;

## convert blast to gff
print "converting blast output $blast_output to gff file\n" if $config{verbose};
my $gff_file = blastGFFConverter($blast_output);

## make bed file from BUSCO file
my $BUSCObed = File::Temp->new( UNLINK => 1, SUFFIX => '.blast' ); # temporary file with BUSCO bed file
open (BUSCOBED, ">$BUSCObed") or die "cannot open file $BUSCObed for writting BUSCO bed file\n";
print "creating Bed from from BUSCO file $config{busco}\n" if $config{verbose};
my %fasta = genometohash($config{busco});
foreach my $seq (keys %fasta) {
	my $seqlen = length $fasta{$seq};
	print BUSCOBED "$seq\t1\t$seqlen\tBUSCO\n";
}
close BUSCOBED;

## calculate coverage
my $coverageBEDoutput = File::Temp->new( UNLINK => 1, SUFFIX => '.bed' ); # temporary file with coverageBED ouptut
print "calculating coverage\n" if $config{verbose};
$cmdoutput = `coverageBed -a $gff_file -b $BUSCObed > $coverageBEDoutput`;
my $totalsize; # size of all BUSCOs combined 
my $coveragesize; # bp covered
open (INPUT, $coverageBEDoutput) or die "cannot open file $coverageBEDoutput\n";
while (my $line = <INPUT>) {
	my @data = split(" ", $line);
	$coveragesize += $data[5];
	$totalsize += $data[6];
}
close INPUT; 

## report output
my $prop = $coveragesize/$totalsize;
my $numberinputseqs = `grep -c ">" $config{fasta}`;
chomp $numberinputseqs;
print "$numberinputseqs\t$prop\n";



###########################################################################
################################ Functions ################################
###########################################################################
sub printError{
    my $string = shift;
    my $exit = shift;
    
    print STDERR $string;
    exit if $exit;
}

###########################################################################
sub printUsage{

print STDOUT "USAGE : RE_estimate_size.pl -f \"reads in fasta format\" -b \"busco file\" -o \"output file\"
Options : 
    -f | fasta		Input fasta file (Mandatory)
    -b | busco		Name of busco file to use (default \"/home/parensburge/db/BUSCO/BUSCO_Vertebrata.fas\"
    -t | threads	number of threads to use (default 12)
    -o | out   		Name of ouput file (default \"out\")
    -v | verbose	MORE text dude !!!!
    -h | help		You already know it ...\n\n";
    exit;
}
