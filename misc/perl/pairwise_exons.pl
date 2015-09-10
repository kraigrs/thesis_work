#!/usr/bin/perl

###################################################################################
#
# 06/07/2010
#
# pairwise_exons.pl
#
# Purpose: from two lists of exons in other species, which ones are common? 
# 
# Input: two lists of exons in other species' space with corresponding dm3 names
#
# Output: list of exons to exclude that are not present in both species 
#
# Usage: perl pairwise_exons.pl <constExonsS1> <constExonsS2> <output>  
#
#        <constExonsS1> ==> a file containing a list of constitutive exons from species 1
#        <constExonsS2> ==> a file containing a list of constitutive exons from species 2
#        <output>       ==> list of exons to exclude
#
###################################################################################

use strict;
use warnings;

main ();
sub main 
{
  my$const = $ARGV[0];
  my$constExonsS1 = $ARGV[1]; 
  my$constExonsS2 = $ARGV[2];
  my$output = $ARGV[3];

  my$exon;
  my@elements;
  my%exons1; my%exons2; my%all_exons;

  open(EXONS,"$const") or die "Can't open $const for reading!\n";
  while(<EXONS>) 
  {
    chomp;
    unless (/^#track.+/)
    {
      @elements = split(/\s+/,$_); 
      $all_exons{$elements[3]} = [$elements[0],$elements[1],$elements[2]];
    }
  }
  close EXONS;

  open(EXONS1,"$constExonsS1") or die "Can't open $constExonsS1 for reading!\n";
  while(<EXONS1>) 
  {
    chomp;
    unless (/^#track.+/)
    {
      @elements = split(/\s+/,$_); 
      $exons1{$elements[3]} = 1;
    }
  }
  close EXONS1;

  open(EXONS2,"$constExonsS2") or die "Can't open $constExonsS2 for reading!\n";
  while(<EXONS2>) 
  {
    chomp;
    unless (/^#track.+/)
    {
      @elements = split(/\s+/,$_); 
      $exons2{$elements[3]} = 1;
    }
  }
  close EXONS2;

  open(OUT,">$output") or die "Error writing to $output!\n";

  foreach $exon (keys %all_exons) 
  {
    unless($exons1{$exon} && $exons2{$exon})
    {
      print OUT "$all_exons{$exon}[0]\t$all_exons{$exon}[1]\t$all_exons{$exon}[2]\t$exon\n";
    }
  }

  close OUT;
}
