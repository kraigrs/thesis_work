#!/usr/bin/perl

########################################################################
# 
# 09/23/2011
#
# line_specific_base_composition.pl
#
# Purpose: for each line in DGRP, number of specific calls
# 
# Input: chromosome- and base-specific line SNP calls
#
# Output: for each line, the number of each base
#
# Syntax: perl compile_SNPs.pl ../DGRP/chr2L_snps.csv ../DGRP/chr2L_SNPs.txt 
#
########################################################################

use strict;
use warnings;

main();

sub main
{
  my$SNP_lines = $ARGV[0];
  my$out = $ARGV[1];

  my$pos; my$call; my$chr; my$line;
  my@elements;
  my%SNPlines; my%calls;

  # compile SNPs from each line to determine a composite value
  print "\nReading in SNP data...\n";
  open(SNP,"$SNP_lines") or die "\nError opening $SNP_lines\n";
  while(<SNP>)
  {
    chomp;
    if(/^chrs/){next;}  
    @elements = split(/,/,$_);
    $chr = $elements[0];
    $pos = $elements[1];
    $line = $elements[2];
    $call = $elements[3];
    $calls{$call} = 1;
    $SNPlines{$line}{$call} += 1;

    #if(!$SNPlines{$line}{$call}){$SNPlines{$line}{$call} = 1;}
    #else{$SNPlines{$line}{$call} += 1;}
  }
  close SNP;
  print "\nDone reading in SNP data, now compiling SNPs...\n";

  open(OUT,">$out") or die "Error writing to $out\n";

  print OUT "line";
  foreach $call (sort keys %calls){print OUT "\t$call";}
 
  foreach $line (keys %SNPlines)
  {
    print OUT "\n$line";
    foreach $call (sort keys %calls)
    {
      if(!$SNPlines{$line}{$call}){print OUT "\t0";}
      else{print OUT "\t$SNPlines{$line}{$call}";}
    }
  }
  close OUT;
}
