#!/usr/bin/perl

########################################################################
# 
# 02/19/2009
#
# changeRef.pl
#
# Purpose: change a reference for later reference by adding
#          a D* suffix (or prefix) i.e. Dmel or Dsec or Dsim
# 
# Input: the original reference
#
# Output: the newly edited reference
#
########################################################################

use strict;
use warnings;

my$fasta = $ARGV[0];
my$output = $ARGV[1];
my$species = $ARGV[2];

open(FASTA,"$fasta") or die "\nError opening $fasta\n";

while(<FASTA>)
{
  chomp;
  if(/^(\>\w+)(\s+.+)$/)
  {
    open(OUT,">> $output") or die "Error writing to $output\n";
    print OUT "$1/$species$2\n";
    close OUT;  
  }
  else
  {
    open(OUT,">> $output") or die "Error writing to $output\n";
    print OUT "$_\n";
    close OUT;
  }
}
close FASTA;
