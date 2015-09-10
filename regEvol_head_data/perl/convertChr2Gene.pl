#!/usr/bin/perl

###################################################################################
#
# 05/16/2010
#
# convertChr2Gene.0.5.pl
#
# Purpose: query a gap file to see if mapped read falls in a gap and, if not, 
#          in which gene and exon it happens to fall in 
# 
# Input: .map file, list of gaps, list of constitutive exons
#
# Output: for every read, an indication of being in a gap, constitutive exon, or neither 
#
# Fix: using intersectBed to MUCH more efficiently intersect BED file with gaps and exons
#
# Usage: perl convertChr2Gene.0.5.pl <bedfile> <gaps> <constExons>   
#        <bedfile>    ==> the .bed file to be converted 
#        <gaps>       ==> a file containing the gap regions
#        <constExons> ==> a file containing a list of constitutive exons
#
###################################################################################

use strict;
use warnings;

main ();
sub main 
{ 
  my$bedfile = $ARGV[0];  
  my$gaps_file = $ARGV[1];
  my$constEx_file = $ARGV[2]; 

  my$cmd = "intersectBed -a $bedfile -b $gaps_file -v | ".
           "intersectBed -a stdin -b $constEx_file -wa -wb -f 1 | ".
           "awk '{print \$10 \"\t\" \$8 \"\t\" \$9 \"\t\" \$4}' ".
           "> $bedfile.converted";

  #print $cmd;

  system $cmd;
}
