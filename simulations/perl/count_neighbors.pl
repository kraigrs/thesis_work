#!/usr/bin/perl

######################################################################################
# 
# 04/27/2013
#
# count_neighbors.pl
#
# Purpose: count the number of "relevant" neighboring differentiating sites
# 
# Input: list of SNPs with positions         
#
# Output: relevant neighbors 
#
######################################################################################

use strict;
use warnings;
use List::Util qw( min max );

my$SNPs = $ARGV[0];
my$overlaps = $ARGV[1];
my$output = $ARGV[2];

my$chr; my$pos; my$n; my$maximum;
my@elements;
my%SNPhash; my%neighbors;

open(SNPS,"$SNPs") or die "\nError opening $SNPs\n";
while(<SNPS>)
{
  chomp;
  @elements = split("\t",$_);

  # mel_mel
  $chr = $elements[0];
  $pos = $elements[2]; # 1-based position in SNP BED file since 1-based in pileup
  $SNPhash{$chr}{$pos} = 1;


  # mel_sim
  #$chr = $elements[0];
  #$pos = $elements[1]; # 1-based position in SNP BED file since 1-based in pileup

}
close SNPS;

open(OVER,"$overlaps") or die "\nError opening $overlaps\n";
while(<OVER>)
{
  chomp;
  @elements = split("\t",$_);

  $chr = $elements[0];
  $n = $elements[3];
  $pos = $elements[6];
  #print "$chr\t$n\t$pos\n"; exit;

  push @{ $neighbors{$chr}{$pos} }, $n;
}
close OVER;

open(OUT,"> $output") or die "Can't open $output to print\n";
print OUT "chr\tpos\trelevant\n";

foreach $chr (keys %SNPhash)
{
  foreach $pos (keys %{ $SNPhash{$chr} })
  {
    if($neighbors{$chr}{$pos})
    {
      $maximum = max @{ $neighbors{$chr}{$pos} };
      print OUT "$chr\t$pos\t$maximum\n";
    }
  }
}
close OUT;
