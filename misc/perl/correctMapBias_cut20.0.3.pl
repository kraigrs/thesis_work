#!/usr/bin/perl

###################################################################################
#
# 02/22/2011
#
# correctMapBias_cut20.0.3.pl
#
# Fix: this file inputs a mixed parental mRNA (cDNA) pool as the parental counts
#
# Purpose: combine summary data from parents and hybrids to obtain common gene set and perform transformation 
# 
# Input: 2 files, 1 summary from mixed parents and 1 summary from hybrid cross
#
# Output: a list of common genes between sets that have been adjusted to species-specific sequencing depth and mapping enrichment 
#
# E.g. perl correctMapBias_cut20.0.3.pl <species1> <species2> <parfile> <hybfile> <out_dir>
#
#      <s1>       ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <s2>       ==> indicate name of 2nd species aligned to
#      <parfile> ==> file with parent counts
#      <hybfile  ==> file with s1 X s2 hybrid counts
#      <out_dir>   ==> directory to send output
#      <suffix>
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
  my$hybfile = $ARGV[3];
  my$out_dir = $ARGV[4];
  my$suffix = $ARGV[5];

  my$line; my$locus; my$s1ct; my$s2ct; my$both; 
  my$parS1; my$parS2;
  my$hybS1; my$hybS2;
  my$factorParS1; my$factorParS2;
  my$factorHybS1; my$factorHybS2;
  my@elements;
  my%loci; my%par; my%hyb;

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
      if($s1ct + $s2ct > 20)
      {
        $loci{$locus} = 1;
        $par{$locus} = [$s1ct,$s2ct,$both]; 
      }
    }
  }
  close FILE;

  # hybfile
  open(FILE,"$hybfile") or die "Can't open $hybfile for reading!";
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
        $hyb{$locus} = [$s1ct,$s2ct,$both];
      } 
    }
  }
  close FILE;

  # Determine accurate counts for each case using common gene sets

#  foreach $locus (keys %loci)   
#  {
#    if($par{$locus} && $hyb{$locus})
#    {  
#      $parS1  += $par{$locus}[0]; $parS2  += $par{$locus}[1];
#      $hybS1 += $hyb{$locus}[0]; $hybS2 += $hyb{$locus}[1];
#    }
#  }

  foreach $locus (keys %loci)   
  {
    if($par{$locus})
    {  
      $parS1  += $par{$locus}[0]; $parS2  += $par{$locus}[1];
    }

    if($hyb{$locus})
    {
      $hybS1 += $hyb{$locus}[0]; $hybS2 += $hyb{$locus}[1];
    }
  }

  # Determine multiplicative factors for sequencing depth and mapping bias adjustment 
  
  # sequencing depth and mapping bias adjustment 
  if($parS1 == $parS2){$factorParS1 = 1; $factorParS2 = 1;}
  elsif($parS1 > $parS2){$factorParS1 = $parS2/$parS1; $factorParS2 = 1;}
  elsif($parS1 < $parS2){$factorParS1 = 1; $factorParS2 = $parS1/$parS2;}
 
  if($hybS1 == $hybS2){$factorHybS1 = 1; $factorHybS2 = 1;}
  elsif($hybS1 > $hybS2){$factorHybS1 = $hybS2/$hybS1; $factorHybS2 = 1;}
  elsif($hybS1 < $hybS2){$factorHybS1 = 1; $factorHybS2 = $hybS1/$hybS2;}

  print "\nfactorParS1 = $factorParS1\nfactorParS2 = $factorParS2\n";
  print "factorHyb1S1 = $factorHybS1\nfactorHybS2 = $factorHybS2\n\n";

  open(PARHYB,">$out_dir/$species1\_$species2.parHyb_mapBias_cut20.$suffix") or die "Error writing to $out_dir/$species1\_$species2.parHyb.$suffix\n";

  print PARHYB  "locus\tPar_$species1\tPar_$species2\tPar_both\tHyb_$species1\tHyb_$species2\tHyb_both";

  foreach $locus (keys %loci)   
  {
    if($par{$locus} && $hyb{$locus})
    {
      printf PARHYB ("\n$locus\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
                      floor($par{$locus}[0]*$factorParS1),floor($par{$locus}[1]*$factorParS2),
                      floor($par{$locus}[2]*$factorParS1+$par{$locus}[2]*$factorParS2),
                      floor($hyb{$locus}[0]*$factorHybS1),floor($hyb{$locus}[1]*$factorHybS2),
                      floor($hyb{$locus}[2]*$factorHybS1+$hyb{$locus}[2]*$factorHybS2));
    }
  }
  close PARHYB;
}
