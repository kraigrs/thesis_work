#!/usr/bin/perl

###################################################################################
#
# 02/23/2010
#
# correctSeqDepth_cut20.0.3.pl
#
# Purpose: combine summary data from parents and hybrids to obtain common gene set and perform transformation 
# 
# Input: 2 files, 1 summary from mixed parents and 1 summary hybrid cross
#
# Output: a list of common genes between sets that have been adjusted to species-specific sequencing depth and mapping enrichment 
#
# E.g. perl correctSeqDepth_cut20.0.2.pl <species1> <species2> <speciesfile> <hyb1file> <hyb1depth> <hyb2file> <hyb2depth> <out_dir>
#
#      <s1>       ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <s2>       ==> indicate name of 2nd species aligned to
#      <parfile> ==> file with s1 parent counts
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
  my$factorS1; my$factorS2; 
  my$factorHyb;
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

  open(PARHYB,">$out_dir/$species1\_$species2.parHyb_seqDepth_cut20.$suffix") or die "Error writing to $out_dir/$species1\_$species2.parHyb.$suffix\n";

  print PARHYB  "locus\tPar_$species1\tPar_$species2\tPar_both\tHyb_$species1\tHyb_$species2\tHyb_both";

  foreach $locus (keys %loci)   
  {
    if($par{$locus} && $hyb{$locus})
    {
      printf PARHYB ("\n$locus\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f",
                      floor($par{$locus}[0]),floor($par{$locus}[1]),floor($par{$locus}[2]),
                      floor($hyb{$locus}[0]),floor($hyb{$locus}[1]),floor($hyb{$locus}[2]));
    }
  }
  close PARHYB;
}
