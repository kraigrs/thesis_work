#!/usr/bin/perl

###################################################################################
#
# 02/12/2010
#
# mergeSAM.pl
#
# Purpose: merge 2 .sam files (ideally from two mates from the same experiment)
# 
# Input: the two .sam files (mate 1 and mate 2). 
#
# Output: a concatenated .sam file containing all reads from the mate pairs (if applicable) 
#
# Example: perl mergeSAM.pl ../mel_sec_data/Hyb_mate1.dmel-gene.sam ../mel_sec_data/Hyb_mate2.dmel-gene.sam ../mel_sec_data/Hyb.dmel-gene.sam &
#
###################################################################################

use strict;
use warnings;

my$SAM_file1 = $ARGV[0];
my$SAM_file2 = $ARGV[1];
my$SAMout = $ARGV[2]; 

open(OUT,"> $SAMout") or die "Error writing to $SAMout\n";

#print "\n\nParsing first mates file...\n";
open(SAM1,"$SAM_file1") or die "\nError opening $SAM_file1\n";
while(<SAM1>)
{
  chomp $_;  
  # make the header
  if(/^@/)
  {
    print OUT "$_\n";
  }
  elsif(/^HWI.+\/1\s+/)
  {
    print OUT "$_\n";
  }
}
close SAM1;

open(SAM2,"$SAM_file2") or die "\nError opening $SAM_file2\n";
while(<SAM2>)
{
  chomp $_;  
  if(/^HWI.+\/2\s+/)
  {
    print OUT "$_\n";
  }
}
close SAM2;
close OUT;
