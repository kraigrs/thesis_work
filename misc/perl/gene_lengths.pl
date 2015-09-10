#!/usr/bin/perl

###################################################################################
#
# 04/17/2012
#
# gene_lengths.pl
#
# Purpose: compile a list of gene names with respective lengths (sum of exon lengths within gene) 
# 
# Input: list of exons (BED format)
#
# Output: list of genes and their lengths 
#
###################################################################################

use strict;
use warnings;

my$const = $ARGV[0];
my$out = $ARGV[1];

our($gene,$start,$stop,$length);
our(@elements);
our(%geneHash);

open(CONST,"$const") or die "Can't open $const for reading!";
while(<CONST>) 
{
  chomp;
  if($_ =~ /^track.+/){next;}
  else
  { 
    @elements = split(/\s/,$_);
    $start = $elements[1];
    $stop = $elements[2];
    $gene = $elements[3];
    $length = $stop-$start;
    $geneHash{$gene} += $length;
  }
}
close CONST;

open(OUT,"> $out") or die "Error writing to $out\n";
print OUT "gene\tlength\n";

foreach $gene (keys %geneHash)
{
  print OUT "$gene\t$geneHash{$gene}\n";
}
close OUT;
