#!/usr/bin/perl
# readdata.pl
# Read in data from both files (orig.data.txt and ill.data.txt) via 
# perl script one line at a time (readdata.pl)
# 
# Store output as single text file, which has TSS and CHR (input.fa)
# Use this file as input to springResearch.pl. springResearch.pl will
# Use downloaded human genome files to output sequence of length 20kb 
# (10kb up and 10kb downstream of TSS) as results.fa

# Read in data from two .txt files 
# and output single text file, which has both TSS# and CHR# for each 
# probe of interest 

use strict;
use warnings;
use springResearchLib;

	# open file for reading
#open(FILE, "orig.data.edited.txt") or die "Could not open file: $!";	
my $mtch = 0;
	# open file for writting
open (FILE3, ">>input.fa") or die ("Cannot open file input.fa \n");
#read probeIDs from orig.data.txt one at a time
	open (FILE, 'orig.data.txt');
	my @orig = <FILE>;
	my @new = 0;
	
	foreach (@orig){
		@new = split('', $_);	
	
		my $probe = $new[2].$new[3].$new[4].$new[5].$new[6].$new[7].$new[8].$new[9];
		#print "probe\n$probe\n";
		open (FILE2, 'ill.data.txt');
		my @array = <FILE2>;
		#my @values = split(' ', $_);

		foreach (@array){
			my @values = split(' ', $_);
			my $illprobe = $values[0];
			$illprobe =~ s/cg//;
			#print "$illprobe\n\n";
			if ($values[0] =~ m/$probe/){
    			printf FILE3 "$_";
    			$mtch = 1;
				last if $mtch ==1;
			}
		
		}

	}
close FILE2;
close FILE3;
close FILE;

