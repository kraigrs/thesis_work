#!/usr/bin/perl

###################################################################################
#
# 03/29/2012
#
# reads_from_BED.pl
#
# Purpose: using read names from a BED file, extract the appropriate reads from FASTA
# 
# Input: BED file, FASTA file 
#
# Output: all reads with matching BED info 
#
# Usage: perl reads_from_BED.pl <BED file> <FASTA file> <output>
#        <BED file>   ==> reads from BED file 
#        <FASTA file>   ==> list of FASTA sequences to extract
#        <output>     ==> output of the new FASTA file
#
###################################################################################

use strict;
use warnings;

my$BED = $ARGV[0];
my$FASTA = $ARGV[1];
my$output = $ARGV[2];

our($name,$seq,$symbol);
our(@elements);
our(%reads);

open(BED,"$BED") or die "\nError opening $BED\n";
while(<BED>)
{
  chomp;  
  @elements = split(/\s+/,$_);
  $name = $elements[3];
  $reads{$name} = 1;
}
close BED;

open(OUT,">$output") or die "Error writing to $output\n";

open(FASTA,"$FASTA") or die "\nError opening $FASTA\n";
while(<FASTA>)
{
  chomp;
  if(/^(\>|\@|\+)(\S+)/){$symbol = $1; $name = $2;} 
  else
  {
    $seq = $_;
    if($reads{$name}){print OUT "$symbol$name\n$seq\n";}
  } 
}
close FASTA;
close OUT;
