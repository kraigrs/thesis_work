#!/usr/bin/perl

###################################################################################
#
# 10/01/2010
#
# correctMapBias_cut20.0.2.pl
#
# Fix: this file inputs a mixed parental mRNA (cDNA) pool as the parental counts
#
# Purpose: combine summary data from parents and hybrids to obtain common gene set and perform transformation 
# 
# Input: 3 files, 1 summary from mixed parents and 1 summary from each reciprocal cross
#
# Output: a list of common genes between sets that have been adjusted to species-specific sequencing depth and mapping enrichment 
#
# E.g. perl processData.pl <species1> <species2> <parfile> <hyb1file> <hyb2file> <out_dir> <suffix>
#
#      <s1>       ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <s2>       ==> indicate name of 2nd species aligned to
#      <parfile> ==> file with mixed parent counts
#      <hyb1file  ==> file with s1 X s2 hybrid counts
#      <hyb2file> ==> file with s2 X s1 hybrid counts
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
  my$parS1; my$parS2;
  my$hyb1S1; my$hyb1S2;
  my$hyb2S1; my$hyb2S2;
  my$factorParS1; my$factorParS2;
  my$factorHyb1S1; my$factorHyb1S2;
  my$factorHyb2S1; my$factorHyb2S2;
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

  # Determine accurate counts for each case using common gene sets

#  foreach $locus (keys %loci)   
#  {
#    if($par{$locus} && $hyb1{$locus} && $hyb2{$locus})
#    {  
#      $parS1  += $par{$locus}[0]; $parS2  += $par{$locus}[1];
#      $hyb1S1 += $hyb1{$locus}[0]; $hyb1S2 += $hyb1{$locus}[1];
#      $hyb2S1 += $hyb2{$locus}[0]; $hyb2S2 += $hyb2{$locus}[1];
#    }
#  }

  foreach $locus (keys %loci)   
  {
    if($par{$locus})
    {  
      $parS1  += $par{$locus}[0]; $parS2  += $par{$locus}[1];
    }

    if($hyb1{$locus} && $hyb2{$locus})
    {
      $hyb1S1 += $hyb1{$locus}[0]; $hyb1S2 += $hyb1{$locus}[1];
      $hyb2S1 += $hyb2{$locus}[0]; $hyb2S2 += $hyb2{$locus}[1];
    }
  }

  # Determine multiplicative factors for sequencing depth and mapping bias adjustment 
  
  # sequencing depth and mapping bias adjustment 
  if($parS1 == $parS2){$factorParS1 = 1; $factorParS2 = 1;}
  elsif($parS1 > $parS2){$factorParS1 = $parS2/$parS1; $factorParS2 = 1;}
  elsif($parS1 < $parS2){$factorParS1 = 1; $factorParS2 = $parS1/$parS2;}
 
  if($hyb1S1 == $hyb1S2){$factorHyb1S1 = 1; $factorHyb1S2 = 1;}
  elsif($hyb1S1 > $hyb1S2){$factorHyb1S1 = $hyb1S2/$hyb1S1; $factorHyb1S2 = 1;}
  elsif($hyb1S1 < $hyb1S2){$factorHyb1S1 = 1; $factorHyb1S2 = $hyb1S1/$hyb1S2;}

  if($hyb2S1 == $hyb2S2){$factorHyb2S1 = 1; $factorHyb2S2 = 1;}
  elsif($hyb2S1 > $hyb2S2){$factorHyb2S1 = $hyb2S2/$hyb2S1; $factorHyb2S2 = 1;}
  elsif($hyb2S1 < $hyb2S2){$factorHyb2S1 = 1; $factorHyb2S2 = $hyb2S1/$hyb2S2;}

  print "\nfactorParS1 = $factorParS1\nfactorParS2 = $factorParS2\n";
  print "factorHyb1S1 = $factorHyb1S1\nfactorHyb1S2 = $factorHyb1S2\n";
  print "factorHyb2S1 = $factorHyb2S1\nfactorHyb2S2 = $factorHyb2S2\n\n";

  open(PARHYB1,">$out_dir/$species1\_$species2.parHyb1_mapBias_cut20.$suffix") or die "Error writing to $out_dir/$species1\_$species2.parHyb1.$suffix\n";
  open(PARHYB2,">$out_dir/$species1\_$species2.parHyb2_mapBias_cut20.$suffix") or die "Error writing to $out_dir/$species1\_$species2.parHyb2.$suffix\n";
  open(HYB1HYB2,">$out_dir/$species1\_$species2.hyb1hyb2_mapBias_cut20.$suffix") or die "Error writing to $out_dir/$species1\_$species2.hyb1hyb2.$suffix\n";

  print PARHYB1  "locus\tPar_$species1\tPar_$species2\tPar_both\tHyb1_$species1\tHyb1_$species2\tHyb1_both";
  print PARHYB2  "locus\tPar_$species1\tPar_$species2\tPar_both\tHyb2_$species1\tHyb2_$species2\tHyb2_both";
  print HYB1HYB2 "locus\tHyb1_$species1\tHyb1_$species2\tHyb1_both\tHyb2_$species1\tHyb2_$species2\tHyb2_both";

  foreach $locus (keys %loci)   
  {
    if($par{$locus} && $hyb1{$locus})
    {
      printf PARHYB1 ("\n$locus\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
                      floor($par{$locus}[0]*$factorParS1),floor($par{$locus}[1]*$factorParS2),
                      floor($par{$locus}[2]*$factorParS1+$par{$locus}[2]*$factorParS2),
                      floor($hyb1{$locus}[0]*$factorHyb1S1),floor($hyb1{$locus}[1]*$factorHyb1S2),
                      floor($hyb1{$locus}[2]*$factorHyb1S1+$hyb1{$locus}[2]*$factorHyb1S2));
    }
    if($par{$locus} && $hyb2{$locus})
    {
      printf PARHYB2 ("\n$locus\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
                      floor($par{$locus}[0]*$factorParS1),floor($par{$locus}[1]*$factorParS2),
                      floor($par{$locus}[2]*$factorParS1+$par{$locus}[2]*$factorParS2),
                      floor($hyb2{$locus}[0]*$factorHyb2S1),floor($hyb2{$locus}[1]*$factorHyb2S2),
                      floor($hyb2{$locus}[2]*$factorHyb2S1+$hyb2{$locus}[2]*$factorHyb2S2));
    }
    if($hyb1{$locus} && $hyb2{$locus})
    {
      printf HYB1HYB2 ("\n$locus\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
                      floor($hyb1{$locus}[0]*$factorHyb1S1),floor($hyb1{$locus}[1]*$factorHyb1S2),
                      floor($hyb1{$locus}[2]*$factorHyb1S1+$hyb1{$locus}[2]*$factorHyb1S2),                      
                      floor($hyb2{$locus}[0]*$factorHyb2S1),floor($hyb2{$locus}[1]*$factorHyb2S2),
                      floor($hyb2{$locus}[2]*$factorHyb2S1+$hyb2{$locus}[2]*$factorHyb2S2));  
    }
  }
  close PARHYB1; close PARHYB2; close HYB1HYB2;
}
