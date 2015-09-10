#!/usr/bin/perl

###################################################################################
#
# 05/23/2010
#
# genomeToExon.pl
#
# Purpose: use IntersectBed to find overlaps to exons from genomic coords 
# 
# Input: bed file of alignments and list of exons (also bed format)
#
# Output: every read that intersects an exon 
#
# Fix: using intersectBed to MUCH more efficiently intersect BED file with exons
#
# Usage: perl genomeToExon.pl <bedfile> <constExons> 
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
           "awk '{print \$8 \"\t\" \$6 \"\t\" \$7 \"\t\" \$4}' ".
           "> $bedfile.converted";

  #print $cmd;
  system $cmd;
}
