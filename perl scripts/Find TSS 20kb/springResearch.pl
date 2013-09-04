#!/usr/bin/perl

#############################################################################
# Brock Donovan
# springResearch.pl
# This perl script finds the sequence 10kb upstream and 10kb downstream 
# for every TSS in the file specified as INPUT. Outputs the sequence as
# results.fa. readdata.pl should be run before this .pl file is executed.
#
# The Illumina reference file used was version hg19
# Reference genome sequence downloaded from UCSC: version hg19 to match up 
# with the Illumina version 36:
# http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/bigZips/

############################################################################# 

#############################################################################

use strict;
use warnings;
#use springResearchLib;

my $TSS;
my $chrom; 
# $kb is the number of bases to retrieve up and down stream
my $kb = 10000;
# open file for reading
open(INPUT, "input.fa") or die "Could not open file: $!";
while (<INPUT>) 
	{
	my @values = split(' ', $_);
    $TSS = $values[2];	
    $chrom = $values[1];
    print "TSS = $TSS";
    print "chrom = $chrom\n";
  	seqfind($TSS, $chrom, $kb);
  	}
close INPUT;



#################################################################################
# seqfind
# Subroutine to find desired nucleotide sequence in human genome
# inputs are $TSS (transcription start site), $chr (which chromosome TSS is on), 
# $kb (how many kb upstream and down stream to retrieve) 
#################################################################################
sub seqfind{
	#my($chrom, $TSS);
	my $TSS = $_[0]; 
	my $chrom = $_[1]; 
	my $kb = $_[2]; 
	my $count = 0;
	
	# open file for reading
	open(FILE, "chr$chrom.fa") or die "Could not open file: $!";
	# open file for writing in append mode >>
	open (FILE2, ">>results.fa") or die ("Cannot open file results.fa \n");
	printf FILE2 "\n>chr$chrom TSS$TSS\n ";
	$count = 0;
	while (<FILE>) 
	{
		my @values = split('', $_);
	
		foreach my $val (@values) {
    		$count = $count +1;
    		#add 5 to $TSS for header characters '>chr_'
    		#grab sequence 10 before and 10 after TSS
    		if ($TSS-$kb <= $count && $count <= $TSS+$kb) {
    		#print "$val";
    		#print "$count\n";
    		printf FILE2 "$val";
    		}      	 
  		}
  		if ($count > $TSS+$kb) {
  				#print "\n";
  				last;
  			}
  		#printf FILE2 "\n\n";
	}
	close FILE;
}

