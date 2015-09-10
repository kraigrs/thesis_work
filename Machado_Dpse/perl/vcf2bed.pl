#!/usr/bin/perl

######################################################################################
# 
# 12/09/2013
#
# vcf2bed.pl
#
# Purpose: read in pileup from SAMtools to call SNPs (or mutations)
# 
# Input: SAMtools alignment pileup and VCF file         
#
# Output: based on certain criteria, return highest confidence variants
#
######################################################################################

use strict;
use warnings;

my$VCF = $ARGV[0];

my$locus; my$pos; my$start;
my@elements;

open(VCF,"$VCF") or die "\nError opening $VCF\n";
while(<VCF>)
{
  chomp;
  unless(/^#/)
  {
    @elements = split("\t",$_);
    $locus = $elements[0];
    $pos = $elements[1];
    $start = $pos-1;
    print "$locus\t$start\t$pos\n";
  }
}
close VCF;
