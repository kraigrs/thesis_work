#!/usr/bin/perl

###################################################################################
#
# 11/07/2012
#
# count_polymorphism.pl
#
# Purpose: count numbers of SNPs in exons and genes 
# 
# Input: list of SNPs by exon
#
# Output: number of SNPs in each exon and gene 
#
###################################################################################

use strict;
use warnings;

my$list = $ARGV[0];
my$exon_out = $ARGV[1];
my$gene_out = $ARGV[2];

my$gene; my$exon; my$start; my$stop; my$length; my$name; my$total; my$SNPs; my$site;
my@elements; my@meta;
my%geneHash;

open(SNP,"$SNPs") or die "Can't open $list for reading!";
while(<SNP>) 
{
  chomp;
  if($_ =~ /^track.+/){next;}
  else
  { 
    @elements = split(/\s/,$_);
    $name = $elements[0];
    @meta = split(/\_/,$name);
    $gene = $meta[0];
    $exon = $meta[1].$meta[2];
    $site = $elements[1];

    $geneHash{$gene}{$exon} += 1;
  }
}
close SNP;

open(OUT1,"> $exon_out") or die "Error writing to $exon_out\n";
open(OUT2,"> $gene_out") or die "Error writing to $gene_out\n";

foreach $gene (keys %geneHash)
{
  $total = 0;

  foreach $exon (keys %{$geneHash{$gene}})
  {
    print OUT1 "$gene\_$exon\t$geneHash{$gene}{$exon}\n";
    $total += $geneHash{$gene}{$exon};
  }

  print OUT2 "$gene\t$total\n";
}
close OUT1;
close OUT2;
