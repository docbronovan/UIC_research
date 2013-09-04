#!/usr/bin/perl
use BioE481Library;
use strict;
use warnings;
#######################################################################################  
# Brock Donovan 10-27-12
# donovan2@uic.edu
# finds a specific protein motif
# uses BioE481Library file
####################################################################################### 
my($protein, $pos, $esiye, $edlya, $esiya, $epiya);
$protein='';
($pos, $esiye, $edlya, $esiya, $epiya) = (0,0,0,0,0);
# extract protein sequence from genbank input gb-sequence.gb
$protein = ExtractProteinFromGenBank();

# motif: E in the 1st position, any in the 2nd position, I or L in 
# the 3rd position, Y in the 4th position, A or E in the 5th position
# Find the position of the matching AA
while ($protein =~ m/E.(I|L)Y(A|E)/g) 
	{
  	print "Found '$&' at character ";
  	#subtract 5 to get the start position of sequence
  	print int(pos($protein))-5;
  	print "\n";
  	#keep count of sequences found
  	if($& eq 'ESIYE') {
  		++$esiye;
  	}
  	if($& eq 'EDLYA') {
  		++$edlya;
  	}
  	if($& eq 'ESIYA') {
  		++$esiya;
  	}
  	if($& eq 'EPIYA'){
  		++$epiya;
  	}
}
print "ESIYE found $esiye times\n";
print "EDLYA found $edlya times\n";
print "ESIYA found $esiya times\n";
print "EPIYA found $epiya times\n";

exit;