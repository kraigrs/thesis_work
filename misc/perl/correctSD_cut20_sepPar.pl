#!/usr/bin/perl

###################################################################################
#
# 03/26/2010
#
# correctSD_cut20_sepPar.pl
#
# Purpose: combine summary data from parents to obtain common gene set and perform transformation 
# 
# Input: one summary file from each parent, aligned to individual respective genome
#
# Output: a list of common genes between sets that have been adjusted to species-specific sequencing depth and mapping enrichment 
#
# E.g. perl correctSD_cut20_sepPar.pl <species1> <species2> <species1file> <s1depth> <species2file> <s2depth> <out_dir>
#
#      <s1>       ==> indicate name of 1st species (i.e. Dmel, Dsec, etc.)
#      <s2>       ==> indicate name of 2nd species
#      <par1file> ==> file with s1 parent counts
#      <s1depth>  ==> number of original sequences in par 1
#      <par2file> ==> file with s2 parent counts
#      <s2depth>  ==> number of original sequences in par2
#      <out_dir>   ==> directory to send output
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
  my$s1depth = $ARGV[3];
  my$par2file = $ARGV[4];
  my$s2depth = $ARGV[5];
  my$out_dir = $ARGV[6];
  my$suffix = $ARGV[7];

  my$line; my$locus; my$s1ct; my$s2ct; 
  my$factorS1; my$factorS2; 
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
      $loci{$locus} = 1;
      $par1{$locus} = $s1ct; 
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
      $s2ct = $elements[1];
      $loci{$locus} = 1;
      $par2{$locus} = $s2ct;
    }
  }
  close FILE;

  # Determine multiplicative factors for sequencing depth and mapping bias adjustment 

  # sequencing depth adjustment
  if($s1depth == $s2depth){$factorS1 = 1; $factorS2 = 1;}
  elsif($s1depth > $s2depth){$factorS1 = $s2depth/$s1depth; $factorS2 = 1;}
  elsif($s1depth < $s2depth){$factorS1 = 1; $factorS2 = $s1depth/$s2depth;}

  print "\nfactorS1 = $factorS1\nfactorS2 = $factorS2\n\n";

  open(PAR,">$out_dir/$species1\_$species2.SD_cut20_sepPar.$suffix") or die "Error writing to $out_dir/$species1\_$species2.SD_cut20_sepPar.$suffix\n";

  print PAR "locus\t$species1\t$species2";

  foreach $locus (keys %loci)   
  {
    if($par1{$locus} && $par2{$locus} && $par1{$locus} + $par2{$locus} >= 20)
    {
      printf PAR ("\n$locus\t%.1f\t%.1f",
                  floor($par1{$locus}*$factorS1),
                  floor($par2{$locus}*$factorS2));
    }
  }
  close PAR;
}
