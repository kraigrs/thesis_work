#!/usr/bin/perl

###################################################################################
#
# 12/06/2011
#
# makeBED_SNPs.pl
#
# Purpose: make BED files for all SNPs identified from DGRP
# 
# Input: SNP files (including the "bad" SNPs) 
#
# Output: SNP files in BED format (to be used later for BEDTools!) 
#
# Usage: perl makeBED_SNPs.pl <SNP file> <bad SNPs> <output>
#        <SNP file>   ==> original SNPs file 
#        <bad SNPs>   ==> list of SNPs that were inconsistent with the reference genome ("bad")
#        <output>     ==> output of the new BED file
#
# E.g. perl makeBED_SNPs.pl ../DGRP/chr2L_SNPs.txt ../DGRP/chr2L_badSNPs.txt ../DGRP/chr2L_SNPs.bed
#
###################################################################################

use strict;
use warnings;

my$SNPs = $ARGV[0];
my$bad = $ARGV[1];
my$out = $ARGV[2];

our($chr,$pos,$code,$posZero);
our(@elements);
our(%SNPhash,%badSNPhash);

open(SNP,"$SNPs") or die "\nError opening $SNPs\n";
while(<SNP>)
{
  chomp;  
  @elements = split(/\s+/,$_);
  $chr = $elements[0];
  $pos = $elements[1];
  $code = $elements[2];
  
  unless($code eq "N"){$SNPhash{$pos} = $code;}
}
close SNP;

open(BAD,"$bad") or die "\nError opening $bad\n";
while(<BAD>)
{
  chomp;  
  @elements = split(/\s+/,$_);
  $chr = $elements[0];
  $pos = $elements[1];
  $code = $elements[2];
  
  $badSNPhash{$pos} = $code;
}
close BAD;

open(OUT,">$out") or die "Error writing to $out\n";
foreach $pos (keys %SNPhash)
{
  unless($badSNPhash{$pos})
  {
    $posZero = $pos - 1;
    print OUT "$chr\t$posZero\t$pos\t$SNPhash{$pos}\n";;
  }
}
close OUT;
