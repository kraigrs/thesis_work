#!/usr/bin/perl

###################################################################################
#
# 03/25/2011
#
# correctMB_cut20_sepPar2.pl
#
# Purpose: combine summary data from parents to obtain common gene set and perform transformation 
# 
# Input: one summary file from each parent, aligned separately to each genome
#
# Output: a list of common genes between sets that have been adjusted to species-specific sequencing depth and mapping enrichment 
#
# E.g. perl correctMB_cut20_sepPar2.pl <species1> <species2> <species1file> <species2file> <out_dir> <suffix>
#
#      <s1>       ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <s2>       ==> indicate name of 2nd species aligned to
#      <par1file> ==> file with s1 parent counts
#      <par2file> ==> file with s2 parent counts
#      <out_dir>  ==> directory to send output
#      <suffix>   ==> what to append on file ("mosaik.exons.txt")
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
  my$par1file = $ARGV[2];
  my$par2file = $ARGV[3];
  my$out_dir = $ARGV[4];
  my$suffix = $ARGV[5];

  my$line; my$locus; my$s1ct; my$s2ct; my$both; 
  my$parS1; my$parS2;
  my$factorParS1; my$factorParS2;

  my@elements; 
  my%loci; my%par1; my%par2;

  # par1file
  open(FILE,"$par1file") or die "Can't open $par1file for reading!";
  LINE:while(<FILE>) 
  {
    chomp;
    $line = $_;  
    if($line =~ /^(gene|exon).+/){next LINE;}
    else
    { 
      @elements = split(/\s/,$line);
      $locus = $elements[0];
      $s1ct = $elements[1];
      $s2ct = $elements[2];
      $both = $elements[3];
      $loci{$locus} = 1;
      $par1{$locus} = [$s1ct,$s2ct,$both]; 
    }
  }
  close FILE;

  # par2file
  open(FILE,"$par2file") or die "Can't open $par2file for reading!";
  LINE:while(<FILE>) 
  {
    chomp;
    $line = $_;  
    if($line =~ /^(gene|exon).+/){next LINE;}
    else
    { 
      @elements = split(/\s/,$line);
      $locus = $elements[0];
      $s1ct = $elements[1];
      $s2ct = $elements[2];
      $both = $elements[3];
      $loci{$locus} = 1;
      $par2{$locus} = [$s1ct,$s2ct,$both];
    }
  }
  close FILE;

  # Determine accurate counts for each case using common gene sets

  foreach $locus (keys %loci)   
  {
    if($par1{$locus} && $par2{$locus} && $par1{$locus}[0] + $par2{$locus}[1] >= 20)
    {
      $parS1  += $par1{$locus}[0]; $parS2  += $par2{$locus}[1];
    }
  }

  # Determine multiplicative factors for sequencing depth and mapping bias adjustment 

  if($parS1 == $parS2){$factorParS1 = 1; $factorParS2 = 1;}
  elsif($parS1 > $parS2){$factorParS1 = $parS2/$parS1; $factorParS2 = 1;}
  elsif($parS1 < $parS2){$factorParS1 = 1; $factorParS2 = $parS1/$parS2;}

  print "\nfactorParS1 = $factorParS1\nfactorParS2 = $factorParS2\n\n";

  open(PAR,">$out_dir/$species1\_$species2.MB_cut20_sepPar2.$suffix") or die "Error writing to $out_dir/$species1\_$species2.MB_cut20_sepPar2.$suffix\n";

  print PAR  "locus\t$species1\t$species2\tboth";

  foreach $locus (keys %loci)   
  {
    if($par1{$locus} && $par2{$locus} && $par1{$locus}[0] + $par2{$locus}[1] >= 20)
    {
      printf PAR ("\n$locus\t%.1f\t%.1f\t%.1f",
                  floor($par1{$locus}[0]*$factorParS1),
                  floor($par2{$locus}[1]*$factorParS2),
                  floor($par1{$locus}[2]*$factorParS1+$par2{$locus}[2]*$factorParS2));
    }   
  }
  close PAR;
}
