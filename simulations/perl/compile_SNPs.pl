#!/usr/bin/perl

########################################################################
# 
# 09/23/2011
#
# compile_SNPs.pl
#
# Purpose: create a consolidated list of SNPs from the DGRP filtered data
# 
# Input: for each chromosome, input the SNPs from each line
#
# Output: no direct out, but saves a condensed list of SNPs by position, basically combining line information
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

  our($pos,$cov,$call,$val,$refArray,$chr,$SNP,$n);
  our(@elements,@alleles);
  our(%SNPlines);

  sub SNP_caller
  {
    my($ref,$bp) = @_;
    my@array = @{$ref};
    my$cons;

    my$A=0;my$C=0;my$G=0;my$T=0;my$N=0;
    my$R=0;my$Y=0;my$M=0;my$K=0;my$S=0;my$W=0;

    foreach my$item (@array)
    {
      if($item eq "A")   {$A = 1;}
      elsif($item eq "C"){$C = 1;}
      elsif($item eq "G"){$G = 1;}
      elsif($item eq "T"){$T = 1;}
      elsif($item eq "N"){$N = 1;}
      elsif($item eq "R"){$R = 1;}
      elsif($item eq "Y"){$Y = 1;}
      elsif($item eq "M"){$M = 1;}
      elsif($item eq "K"){$K = 1;}
      elsif($item eq "S"){$S = 1;}
      elsif($item eq "W"){$W = 1;}
      else{die "Error: SNP call encountered that is not on the list!";}
    }
    my$sumNuc = $A+$C+$G+$T;
    my$sumCode = $R+$Y+$M+$K+$S+$W;

    if($sumNuc > 2 || $sumCode > 1)
    {
      #print "\nLine values at $bp inconsistent with SNP call!";
      #print "\nR=$R\tA=$A\tG=$G";
      #print "\nY=$Y\tC=$C\tT=$T";
      #print "\nM=$M\tA=$A\tC=$C";
      #print "\nK=$K\tG=$G\tT=$T";
      #print "\nS=$S\tC=$C\tG=$G";
      #print "\nW=$W\tA=$A\tT=$T\n";
      $cons = "N";
    }
    elsif($sumCode == 1)
    {
      if($R)
      {
        if($C || $T){$cons = "N";}
        else{$cons = "R";}
      }
      elsif($Y)
      {
        if($A || $G){$cons = "N";}
        else{$cons = "Y";}
      }
      elsif($M)
      {
        if($T || $G){$cons = "N";}
        else{$cons = "M";}
      }
      elsif($K)
      {
        if($A || $C){$cons = "N";}
        else{$cons = "K";}
      }
      elsif($S)
      {
        if($A || $T){$cons = "N";}
        else{$cons = "S";}
      }
      elsif($W)
      {
        if($G || $C){$cons = "N";}
        else{$cons = "W";}
      }
    }
    elsif($sumNuc == 2)
    {
      if($A && $G){$cons = "R";}
      elsif($C && $T){$cons = "Y";}
      elsif($A && $C){$cons = "M";}
      elsif($T && $G){$cons = "K";}
      elsif($C && $G){$cons = "S";}
      elsif($A && $T){$cons = "W";}
    }
    elsif($sumNuc == 1)
    {
      if($A){$cons = "N";}
      elsif($C){$cons = "N";}
      elsif($G){$cons = "N";}
      elsif($T){$cons = "N";}
    }
    else{die "Error: Line values inconsistent with SNP call!";}

    return $cons;
  }

  # compile SNPs from each line to determine a composite value
  print "\nReading in SNP data...\n";
  open(SNP,"$SNP_lines") or die "\nError opening $SNP_lines\n";
  HERE:while(<SNP>)
  {
    chomp;
    if(/^chrs/){next HERE;}  
    @elements = split(/,/,$_);
    $chr = $elements[0];
    $pos = $elements[1];
    $call = $elements[3];

    push @{ $SNPlines{$pos} },$call;
    #$n = scalar(@{ $SNPlines{$pos} }); 
    #print "\n$n";
  }
  close SNP;
  print "\nDone reading in SNP data, now compiling SNPs...\n";

  open(OUT,">$out") or die "Error writing to $out\n";
  foreach $pos (keys %SNPlines)
  {
    $refArray = \@{ $SNPlines{$pos} };
    #print "\nNo error in referencing array\n";
    $SNP = &SNP_caller($refArray,$pos);
    #print "\nNo error in calling subroutine\n";
    unless($SNP eq "N"){print OUT "$chr\t$pos\t$SNP\n";}
  }
  close OUT;
}
