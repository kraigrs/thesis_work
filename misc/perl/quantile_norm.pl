#!/usr/bin/perl

###################################################################################
#
# 06/15/2011
#
# quantile_norm.pl
#
# Purpose: quantile normalize parental and hybrid data 
# 
# Input: 2 files, parents and hybrids or both hybrid directions
#
# Output: quantile normalized counts 
#
# E.g. perl quantile_norm.pl <species1> <species2> <file1> <file2> <out_dir>
#
#      <species1> ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <species2> ==> indicate name of 2nd species aligned to
#      <file1>    ==> file with s1 parent counts
#      <file2>    ==> file with s2 parent counts
#      <out_dir>  ==> directory to send output
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
  my$species1 = $ARGV[0]; 
  my$species2 = $ARGV[1];
  my$file1 = $ARGV[2];
  my$file2 = $ARGV[3];
  my$out_dir = $ARGV[4];

  my$line; my$locus; my$s1ct; my$s2ct;
  my@elements; my@list; my@file1species1; my@file1species2; my@file2species1; my@file2species2;
  my%loci; 

  # Quantile normalization
  # Usage: quantNorm(<file1species1>,<file1species2>,<file2species1>,<file2species2>) = quantile-normalized matrix in original order
  sub quantNorm
  {
    my$arrayRef = @_;
    my@sorted = sort {$a <=> $b} @{$arrayRef};
    my$N = @sorted;
    my$h = ($N-1)*$p + 1;
    my$Qp = $sorted[floor($h)-1] + ($h - floor($h))*($sorted[floor($h)] - $sorted[floor($h)-1]);
    return $Qp;
  }

  # file1
  open(FILE,"$file1") or die "Can't open $file1 for reading!";
  LINE:while(<FILE>) 
  {
    chomp;
    $line = $_;  
    if($line =~ /^(gene|exon).+/){next LINE;}
    elsif($line =~ /^CG.+/)
    { 
      @elements = split(/\s/,$line);
      $locus = $elements[0];
      $s1ct = $elements[1];
      $s2ct = $elements[2];
      
      if($s1ct + $s2ct > 20)
      {
        $loci{$locus} = 1;
        $file1{$locus} = [$s1ct,$s2ct]; 
      }
    }
  }
  close FILE;

  # file2
  open(FILE,"$file2") or die "Can't open $file2 for reading!";
  LINE:while(<FILE>) 
  {
    chomp;
    $line = $_;  
    if($line =~ /^(gene|exon).+/){next LINE;}
    elsif($line =~ /^CG.+/)
    { 
      @elements = split(/\s/,$line);
      $locus = $elements[0];
      $s1ct = $elements[1];
      $s2ct = $elements[2];
      
      if($s1ct + $s2ct > 20)
      {
        $loci{$locus} = 1;
        $file2{$locus} = [$s1ct,$s2ct]; 
      }
    }
  }
  close FILE;

  foreach $locus (keys %loci)   
  {
    if($file1{$locus} && $file2{$locus})
    {
      push @file1species1,$file1{$locus}[0];
      push @file1species2,$file1{$locus}[1];

      push @file2species1,$file2{$locus}[0];
      push @file2species2,$file2{$locus}[1];

      push @list,$locus;
    }
  }

}
