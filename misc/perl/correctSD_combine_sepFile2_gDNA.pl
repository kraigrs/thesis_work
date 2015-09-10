#!/usr/bin/perl

###################################################################################
#
# 01/10/2012
#
# correctSD_cut20_combine_sepFile2.pl
#
# Purpose: correct counts from two sequencing files to adjust for sequencing depth 
# 
# Input: two summary files, aligned to 2 genomes
#
# Output: a list of common genes between the files that have been adjusted to sequencing depth
#
# E.g. perl correctSD_cut20_sepPar2.pl <species1> <species2> <name1> <file1> <file1depth> <name2> <file2> <file2depth> <output>
#
#      <species1> ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <species2> ==> indicate name of 2nd species aligned to
#      <name1> ==> name of file 1 
#      <file1> ==> file with s1, s2, and Both counts
#      <file1depth>  ==> number of original sequences in file1
#      <name2> ==> name of file 2 
#      <file2> ==> file with s1, s2, and Both counts
#      <file2depth>  ==> number of original sequences in file2
#      <output> ==> should be clear (e.g. zhr+z30_zhrXz30_SD_cut20.mosaik.exons.txt)
#
###################################################################################

use strict;
use warnings;
#use POSIX;

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
  my$name1 = $ARGV[2];
  my$par1file = $ARGV[3];
  my$s1depth = $ARGV[4];
  my$name2 = $ARGV[5];
  my$par2file = $ARGV[6];
  my$s2depth = $ARGV[7];
  my$output = $ARGV[8];

  my$line; my$locus; my$s1ct; my$s2ct; my$both;
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

  # Determine multiplicative factors for sequencing depth and mapping bias adjustment 

  # sequencing depth adjustment
  if($s1depth == $s2depth){$factorS1 = 1; $factorS2 = 1;}
  elsif($s1depth > $s2depth){$factorS1 = $s2depth/$s1depth; $factorS2 = 1;}
  elsif($s1depth < $s2depth){$factorS1 = 1; $factorS2 = $s1depth/$s2depth;}

  print "\nfactorS1 = $factorS1\nfactorS2 = $factorS2\n\n";

  open(OUT,">$output") or die "Error writing to $output\n";

  # separate parents
  print OUT "locus\t$name1\_$species1\t$name1\_$species2\t$name1\_Both\t$name2\_$species1\t$name2\_$species2\t$name2\_Both";

  foreach $locus (keys %loci)   
  {
    # use this if combining separate parents

    if($par1{$locus} && $par2{$locus})
    {
      printf OUT ("\n$locus\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f",
                  int($par1{$locus}[0]*$factorS1+0.5),
                  int($par1{$locus}[1]*$factorS1+0.5),
                  int($par1{$locus}[2]*$factorS1+0.5),
                  int($par2{$locus}[0]*$factorS2+0.5),
                  int($par2{$locus}[1]*$factorS2+0.5),
                  int($par2{$locus}[2]*$factorS2+0.5))
    }
  }
  close OUT;
}
