#!/usr/bin/perl

########################################################################
# 
# 07/16/2010
#
# GEO2FASTA.pl
#
# Purpose: makes a FASTA file from data downloaded from GEO
# 
# Input: a file containing read, tile, lane, mate, qual
#
# Output: FASTA file
#
########################################################################

use strict;
use warnings;

main ();
sub main 
{
  my$file = $ARGV[0];
  my$mate = $ARGV[1];
  my$output = $ARGV[2];

  my$name;
  my$found = 0;
  my$seq;

  open(FILE,"$file") or die "\nError opening $file\n";
  open(OUT,">$output") or die "Error writing to $output\n";

  while(<FILE>)
  {
    chomp;
    #print "$_\n";
    if(/^\@\S+\s+(HWI\S+)\s+\S+$/)
    {
      $name = $1;
      $found = 1;
    }
    elsif($found == 1)
    { 
      $seq = $_;    
      print OUT "\>$name\/$mate\n$seq\n";
      $found = 0;
    }
  }
  close FILE;
}
