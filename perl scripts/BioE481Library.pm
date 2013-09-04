#!/usr/bin/perl

############################################### 
# Brock Donovan
# PERL LIBRARY for BioE projects
###############################################

#################################################################################
# BaseCount
# Subroutine to count A's, T's, C's, G's
#################################################################################
sub BaseCount{
	while($dna =~ m/a/ig){$a++};
	while($dna =~ m/t/ig){$t++};
	while($dna =~ m/c/ig){$c++};
	while($dna =~ m/g/ig){$g++};
	while($dna =~ m/n/ig){$n++};

	return ($a,$c,$t,$g,$n);
}
#################################################################################
# CheckBase
# Check input sequence for invalid characters
#################################################################################
sub CheckBase 
{
  my($dna) = @_;
  #Convert the DNA from scalar string form to array form
  my(@dna) = split ('', $dna);
  foreach my $base (@dna) 
  {
    unless($base =~/[ATGCN]/ig)
    {
      print "\n$base is not a valid nucleotide base!\n\n";
      exit;
    }
  }
}
#################################################################################
# CheckProtein
# Subroutine to check amino acid validity
#################################################################################
sub CheckProtein{
	unless($protein =~ m/[AVWIMFLPGSTCYNQDEKRH]/ig){
		print "\n$protein is not a valid protein!\n\n";}
	exit;
}
#################################################################################
# codon2aa
# Creates gene to AA hash
#################################################################################
sub CodonToaa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = 
    (
	'TCA' => 'S',	# Serine
	'TCC' => 'S',	# Serine
	'TCG' => 'S',	# Serine
	'TCT' => 'S',	# Serine
	'TTC' => 'F',	# Phenylalanine
	'TTT' => 'F',	# Phenylalanine
	'TTA' => 'L',	# Leucine
	'TTG' => 'L',	# Leucine
	'TAC' => 'Y',	# Tyrosine
	'TAT' => 'Y',	# Tyrosine
	'TAA' => '_',	# Stop
	'TAG' => '_',	# Stop
	'TGC' => 'C',	# Cysteine
	'TGT' => 'C',	# Cysteine
	'TGA' => '_',	# Stop
	'TGG' => 'W',	# Tryptophan
	'CTA' => 'L',	# Leucine
	'CTC' => 'L',	# Leucine
	'CTG' => 'L',	# Leucine
	'CTT' => 'L',	# Leucine
	'CCA' => 'P',	# Proline
	'CCC' => 'P',	# Proline
	'CCG' => 'P',	# Proline
	'CCT' => 'P',	# Proline
	'CAC' => 'H',	# Histidine
	'CAT' => 'H',	# Histidine
	'CAA' => 'Q',	# Glutamine
	'CAG' => 'Q',	# Glutamine
	'CGA' => 'R',	# Arginine
	'CGC' => 'R',	# Arginine
	'CGG' => 'R',	# Arginine
	'CGT' => 'R',	# Arginine
	'ATA' => 'I',	# Isoleucine
	'ATC' => 'I',	# Isoleucine
	'ATT' => 'I',	# Isoleucine
	'ATG' => 'M',	# Methionine
	'ACA' => 'T',	# Threonine
	'ACC' => 'T',	# Threonine
	'ACG' => 'T',	# Threonine
	'ACT' => 'T',	# Threonine
	'AAC' => 'N',	# Asparagine
	'AAT' => 'N',	# Asparagine
	'AAA' => 'K',	# Lysine
	'AAG' => 'K',	# Lysine
	'AGC' => 'S',	# Serine
	'AGT' => 'S',	# Serine
	'AGA' => 'R',	# Arginine
	'AGG' => 'R',	# Arginine
	'GTA' => 'V',	# Valine
	'GTC' => 'V',	# Valine
	'GTG' => 'V',	# Valine
	'GTT' => 'V',	# Valine
	'GCA' => 'A',	# Alanine
	'GCC' => 'A',	# Alanine
	'GCG' => 'A',	# Alanine
	'GCT' => 'A',	# Alanine
	'GAC' => 'D',	# Aspartic Acid
	'GAT' => 'D',	# Aspartic Acid
	'GAA' => 'E',	# Glutamic Acid
	'GAG' => 'E',	# Glutamic Acid
	'GGA' => 'G',	# Glycine
	'GGC' => 'G',	# Glycine
	'GGG' => 'G',	# Glycine
	'GGT' => 'G',	# Glycine
    );

    if(exists $genetic_code{$codon}) 
    {
    	return $genetic_code{$codon};
    }
    else
    {
    	print STDERR "Bad codon \"$codon\"!!\n";
    	exit;
    }
}
#################################################################################
# ExtractDNAFromGenBank
# subroutine that returns the gene sequence from input file gb-sequence.gb
#################################################################################
sub ExtractDNAFromGenBank{
	my(@array,$line,$dna,$flag);
	$flag=0;
	$dna='';
	# read in file into @array
	open (FILE, 'gb-sequence.gb');
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
	return $dna;
	exit;
}
#################################################################################
# ExtractProteinFromGenBank
# subroutine that returns the gene sequence from input file gb-sequence.gb
#################################################################################
sub ExtractProteinFromGenBank
{
	my(@array,$line,$protein,$flag);
	$flag=0;
	$protein='';
	# read in file into @array
	open (FILE, 'gb-sequence.gb');
	@array = <FILE>;
	#	look through each value in the array for protein sequence start point,
	#	then store result in $protein
	foreach (@array){
	#	set flag =1 if come across /translation="
		#print $flag;
		#print $_;
		if($_ =~ m/\/translation=\"/i){
			#$protein=$_;
			$flag = 1;
		}
	#	set flag = 0 if come across exon	
		if($_ =~ m/exon/){
			$flag=0;
		}	
	#	if flag = 1 concatenate current value from array to $protein
		if($flag==1){
			$protein=$protein.$_;
			#print $protein;
		}
	#	remove translation and anything thats not letters from $protein
		$protein =~ s/[^a-z]//ig;
		$protein =~ s/translation//ig;
	}	
	return $protein;
	exit;
}

#################################################################################
# ReverseComplement
# Subroutine that returns the reverse complement of a sequence
#################################################################################
sub ReverseComplement{
	my$s1 = uc $s1;
	$s1 =~ tr/[A,T,C,G]/[T,A,G,C]/;
		return($s1);
	exit;
}
#################################################################################
# Translate
# Subroutine Translate 
#################################################################################
sub Translate{
	my($dna,$codon,$x,$aa,$remainder); 
	$x=0;
	my($dna) = @_;
	$remainder = length($dna) % 3;
	# print "remainder $remainder\n";
	print "Nascent protein string: ";
	while ($x<(length($dna)-$remainder)){
		$codon = substr($dna, $x, 3);
		$aa = $aa.CodonToaa($codon);
		$x = $x+3;
	}
	return($aa); 
	exit;
}


#################################################################################
# WriteToFile
# Subroutine which is sent a scalar of text, and then query 
# the user for a filename to write the file to. Will open the 
# filename and write the text to it, then close the file. 
#################################################################################
sub WriteToFile{
	print "Please enter a file name to write the data to: ";
	chomp(my$file=<>);
	open FILE, '>$file.txt' or die "Could not open file inputfile: $!";
	close FILE;
exit;
}



1;