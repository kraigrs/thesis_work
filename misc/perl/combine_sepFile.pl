#!/usr/bin/perl

###################################################################################
#
# 01/10/2012
#
# correctSD_combine_sepFile2.pl
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

main();
sub main
{
  my$species1 = $ARGV[0]; 
  my$species2 = $ARGV[1];
  my$par1file = $ARGV[2];
  my$par2file = $ARGV[3];
  my$output = $ARGV[4];

  my$line; my$locus; my$s1ct; my$s2ct;
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

  open(OUT,">$output") or die "Error writing to $output\n";

  # separate parents
  print OUT "locus\t$species1\t$species2";

  foreach $locus (keys %loci)   
  {
    # use this if combining separate parents
    if($par1{$locus} && $par2{$locus})
    {
      print OUT "\n$locus\t$par1{$locus}\t$par2{$locus}";
    }
  }
  close OUT;
}
