#!/usr/bin/perl

###################################################################################
#
# 05/23/2010
#
# convertChr2Gene.0.6.pl
#
# Purpose: query an exon file to see if mapped read intersects it 
# 
# Input: .map file, list of constitutive exons
#
# Output: every read that intersects an exon 
#
# Fix: using intersectBed to MUCH more efficiently intersect BED file with exons
#
# Usage: perl convertChr2Gene.0.6.pl <bedfile> <constExons> 
#        <bedfile>    ==> the .bed file to be converted 
#        <constExons> ==> a file containing a list of constitutive exons
#
###################################################################################

use strict;
use warnings;

main ();
sub main 
{ 
  my$bedfile = $ARGV[0];  
  my$constEx_file = $ARGV[1]; 

  my$cmd = "intersectBed -a $bedfile -b $constEx_file -wa -wb -f 1 | ".
           "awk '{print \$13 \"\t\" \$11 \"\t\" \$12 \"\t\" \$4}' ".
           "> $bedfile.converted";

  #print $cmd;
  system $cmd;
}
