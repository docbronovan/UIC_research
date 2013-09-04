#!/usr/bin/perl

use strict;
use warnings;
#######################################################################################  
# Brock Donovan 2-18-13
# donovan2@uic.edu
# NU project
# 
# read in file, parse through each line, extract DNA sequence
####################################################################################### 
my(@array,$line,$dna,$flag);
$flag=0;
$dna='';
# read in file into @array
open (FILE, 'test.fa');
 while (<FILE>) {
 	print "$_";
 	# split into an array
 }
@array = <FILE>;

#	look through each value in the array for DNA sequence start point,
#	then store result in $dna
foreach (@array){
#	set flag =1 if come across ORIGIN
	if($_ =~ m/ORIGIN/i){
		$flag=1;
		next;
		}
#	set flag = 0 if come across //	
	if($_ =~ m/\/\//){
		$flag=0;
		}	
#	if flag = 1 concatenate current value from array to $dna
	if($flag==1){
		$dna=$dna.$_;
		}
#	remove anything thats not letters from $dna
	$dna =~ s/[^a-z]//ig;
	}	
print $dna;
exit;