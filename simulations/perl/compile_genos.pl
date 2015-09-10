#!/usr/bin/perl

########################################################################
# 
# 09/23/2011
#
# compile_genos.pl
#
# Purpose: for a single line from DGRP, filter out heterozygous sites and N calls
# 
# Input: line SNP info
#
# Output: true DGRP genotype calls
#
# Syntax: perl compile_genos.pl ../DGRP/DGRP_line_40_snps.csv ../DGRP/DGRP_line_40_genos.bed 
#
########################################################################

use strict;
use warnings;

main();

sub main
{
  my$line = $ARGV[0];
  my$out = $ARGV[1];

  my$pos; my$call; my$chr; my$zero;
  my@elements;

  open(OUT,">$out") or die "Error writing to $out\n";

  open(LINE,"$line") or die "\nError opening $line\n";
  while(<LINE>)
  {
    chomp;
    if(/^chrs/){next}  
    @elements = split(/,/,$_);
    $chr = $elements[0];
    $pos = $elements[1];
    $call = $elements[3];
    $zero = $pos - 1;

    if($call eq "A"){print OUT "chr$chr\t$zero\t$pos\t$call\n";}
    elsif($call eq "C"){print OUT "chr$chr\t$zero\t$pos\t$call\n";}
    elsif($call eq "G"){print OUT "chr$chr\t$zero\t$pos\t$call\n";}
    elsif($call eq "T"){print OUT "chr$chr\t$zero\t$pos\t$call\n";}
  }
  close LINE;
  close OUT;
}
