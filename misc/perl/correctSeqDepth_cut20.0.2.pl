#!/usr/bin/perl

###################################################################################
#
# 10/01/2010
#
# correctSeqDepth_cut20.0.2.pl
#
# Purpose: combine summary data from parents and hybrids to obtain common gene set and perform transformation 
# 
# Input: 3 files, 1 summary from mixed parents and 1 summary from each reciprocal cross
#
# Output: a list of common genes between sets that have been adjusted to species-specific sequencing depth and mapping enrichment 
#
# E.g. perl correctSeqDepth_cut20.0.2.pl <species1> <species2> <speciesfile> <hyb1file> <hyb1depth> <hyb2file> <hyb2depth> <out_dir>
#
#      <s1>       ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <s2>       ==> indicate name of 2nd species aligned to
#      <parfile> ==> file with mixed parent counts
#      <hyb1file  ==> file with s1 X s2 hybrid counts
#      <hyb1depth>  ==> number of original sequences in hyb 1
#      <hyb2file> ==> file with s2 X s1 hybrid counts
#      <hyb1depth>  ==> number of original sequences in hyb 1
#      <out_dir>   ==> directory to send output
#      <suffix>  ==> something like "mosaik.txt", anything you want to add to end of file names
#  
###################################################################################

use strict;
use warnings;
use POSIX;

#use Statistics::Descriptive;
#$stat = Statistics::Descriptive::Full->new();
#$stat->add_data(1,2,3,4); $mean = $stat->mean();
#$var  = $stat->variance();
#$tm   = $stat->trimmed_mean(.25);
#$Statistics::Descriptive::Tolerance = 1e-10;

main();
sub main
{
  my$species1 = $ARGV[0]; 
  my$species2 = $ARGV[1];
  my$parfile = $ARGV[2];
  my$hyb1file = $ARGV[3];
  my$hyb1depth = $ARGV[4];
  my$hyb2file = $ARGV[5];
  my$hyb2depth = $ARGV[6];
  my$out_dir = $ARGV[7];
  my$suffix = $ARGV[8];

  my$line; my$locus; my$s1ct; my$s2ct; my$both; 
  my$factorS1; my$factorS2; 
  my$factorHyb1; my$factorHyb2;
  my@elements; my@par1cts; my@par2cts; 
  my%loci; my%par; my%hyb1; my%hyb2;

  # parfile
  open(FILE,"$parfile") or die "Can't open $parfile for reading!";
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
      $both = $elements[3];
      push @par1cts,$s1ct;
      if($s1ct + $s2ct > 20)
      {
        $loci{$locus} = 1;
        $par{$locus} = [$s1ct,$s2ct,$both]; 
      }
    }
  }
  close FILE;

  # hyb1file
  open(FILE,"$hyb1file") or die "Can't open $hyb1file for reading!";
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
      $both = $elements[3];
      if($s1ct + $s2ct > 20)
      {
        $loci{$locus} = 1;
        $hyb1{$locus} = [$s1ct,$s2ct,$both];
      } 
    }
  }
  close FILE;

  # hyb2file
  open(FILE,"$hyb2file") or die "Can't open $hyb2file for reading!";
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
      $both = $elements[3];
      if($s1ct + $s2ct > 20)
      {
        $loci{$locus} = 1;
        $hyb2{$locus} = [$s1ct,$s2ct,$both]; 
      }
    }
  }
  close FILE;

  # Determine multiplicative factors for sequencing depth and mapping bias adjustment 

  # sequencing depth adjustment (hybrids only)
  
  if($hyb1depth == $hyb2depth){$factorHyb1 = 1; $factorHyb2 = 1;}
  elsif($hyb1depth > $hyb2depth){$factorHyb1 = $hyb2depth/$hyb1depth; $factorHyb2 = 1;}
  elsif($hyb1depth < $hyb2depth){$factorHyb1 = 1; $factorHyb2 = $hyb1depth/$hyb2depth;}

  print "\nfactorHyb1 = $factorHyb1\nfactorHyb2 = $factorHyb2\n\n";

  open(PARHYB1,">$out_dir/$species1\_$species2.parHyb1_seqDepth_cut20.$suffix") or die "Error writing to $out_dir/$species1\_$species2.parHyb1.$suffix\n";
  open(PARHYB2,">$out_dir/$species1\_$species2.parHyb2_seqDepth_cut20.$suffix") or die "Error writing to $out_dir/$species1\_$species2.parHyb2.$suffix\n";
  open(HYB1HYB2,">$out_dir/$species1\_$species2.hyb1hyb2_seqDepth_cut20.$suffix") or die "Error writing to $out_dir/$species1\_$species2.hyb1hyb2.$suffix\n";

  print PARHYB1  "locus\tPar_$species1\tPar_$species2\tPar_both\tHyb1_$species1\tHyb1_$species2\tHyb1_both";
  print PARHYB2  "locus\tPar_$species1\tPar_$species2\tPar_both\tHyb2_$species1\tHyb2_$species2\tHyb2_both";
  print HYB1HYB2 "locus\tHyb1_$species1\tHyb1_$species2\tHyb1_both\tHyb2_$species1\tHyb2_$species2\tHyb2_both";

  foreach $locus (keys %loci)   
  {
    if($par{$locus} && $hyb1{$locus})
    {
      printf PARHYB1 ("\n$locus\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
                      floor($par{$locus}[0]),floor($par{$locus}[1]),floor($par{$locus}[2]),
                      floor($hyb1{$locus}[0]),floor($hyb1{$locus}[1]),floor($hyb1{$locus}[2]));
    }
    if($par{$locus} && $hyb2{$locus})
    {
      printf PARHYB2 ("\n$locus\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
                      floor($par{$locus}[0]),floor($par{$locus}[1]),floor($par{$locus}[2]),
                      floor($hyb2{$locus}[0]),floor($hyb2{$locus}[1]),floor($hyb2{$locus}[2]));
    }
    if($hyb1{$locus} && $hyb2{$locus})
    {
      printf HYB1HYB2 ("\n$locus\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
                      floor($hyb1{$locus}[0]*$factorHyb1),floor($hyb1{$locus}[1]*$factorHyb1),floor($hyb1{$locus}[2]*$factorHyb1),
                      floor($hyb2{$locus}[0]*$factorHyb2),floor($hyb2{$locus}[1]*$factorHyb2),floor($hyb2{$locus}[2]*$factorHyb2));
    }
  }
  close PARHYB1; close PARHYB2; close HYB1HYB2;
}
