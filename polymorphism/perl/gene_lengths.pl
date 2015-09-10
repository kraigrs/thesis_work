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
#my$exon_out = $ARGV[1];
#my$gene_out = $ARGV[2];

my$gene; my$exon; my$start; my$stop; my$length; my$name; my$total;
my@elements; my@meta;
my%geneHash;

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
    $exon = $start."_".$stop;
    $name = $elements[3];
    @meta = split(/\_/,$name);
    $gene = $meta[0];
    #$gene = $elements[3];
    $length = $stop-$start;

    $geneHash{$gene}{$exon} += $length;
  }
}
close CONST;

#open(OUT1,"> $exon_out") or die "Error writing to $exon_out\n";
#open(OUT2,"> $gene_out") or die "Error writing to $gene_out\n";

foreach $gene (keys %geneHash)
{
  $total = 0;

  foreach $exon (keys %{$geneHash{$gene}})
  {
    #print OUT1 "$gene\_$exon\t$geneHash{$gene}{$exon}\n";
    $total += $geneHash{$gene}{$exon};
  }

  print "$gene\t$total\n";
}
#print "\n$total\n\n";
#close OUT1;
#close OUT2;
