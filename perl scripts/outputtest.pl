#!/usr/bin/perl
# Brock Donovan 12-2-12
# Test script. Checks how many lines, characters, and words are a file

use strict;
use warnings;

open(FILE, "chr21.fa") or die "Could not open file: $!";

my ($lines, $words, $chars) = (0,0,0);

while (<FILE>) {
    $lines++;
    $chars += length($_);
    $words += scalar(split(/\s+/, $_));
}

print("lines=$lines words=$words chars=$chars\n");
close FILE;

#Create output1.txt file with a new one						
open (FILE, ">output1.txt") or die ("Cannot open file output1.txt \n");
printf FILE "lines=$lines words=$words chars=$chars\n";
close FILE;
