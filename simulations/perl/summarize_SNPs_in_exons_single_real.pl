#!/usr/bin/perl

######################################################################################
# 
# 12/22/2012
#
# summarize_SNPs_in_exons_single.pl
#
# Purpose: for a set of SNPs, how many are in each exon, how many have >= 20 reads
# 
# Input:        
#
# Output: 
#
######################################################################################

use strict;
use warnings;

my$SNPs = $ARGV[0];
my$mm = $ARGV[1];

my$chr; my$pos; my$total; my$density; my$map; my$exon;
my$n; my$k; my$unbiased_n; my$unbiased_k; my$perfect_n; my$perfect_k;
my@elements;
my%set; my%unbiased; my%perfect;

open(SNPS,"$SNPs") or die "\nError opening $SNPs\n";
while(<SNPS>)
{
  chomp;
  @elements = split("\t",$_);
  $exon = $elements[0];
  $pos = $elements[1]; # 1-based position
  $total = $elements[2]; # number of overlapping reads
  $density = $elements[9]; # density of differentiating sites
  $map = $elements[10]/$elements[11]; # mappability

  $set{$exon}{$pos} = [$total,$density,$map];
}
close SNPS;

$n = 0; $k = 0;
$unbiased_k = 0; 
$perfect_k = 0;

foreach $exon (keys %set)
{
  $n += 1;

  foreach $pos (keys %{$set{$exon}})
  {
    $k += 1;

    if($set{$exon}{$pos}[1] <= $mm)
    {
      $unbiased_k += 1;
      $unbiased{$exon} = 1;

      if($set{$exon}{$pos}[2] == 1)
      {
        $perfect_k += 1;
        $perfect{$exon} = 1; 
      }
    }
  }
}

$unbiased_n = keys %unbiased;
$perfect_n = keys %perfect;

# summary

print "\nNumber of SNPs: $k\n"; 
print "Number of exons with > 0 SNPs : $n\n\n"; 
print "\tNumber of unbiased SNPs: $unbiased_k\n";
print "\tNumber of exons with > 0 unbiased SNPs: $unbiased_n\n\n";
print "\t\tNumber of unbiased and perfectly-mappable SNPs: $perfect_k\n";
print "\t\tNumber of exons with > 0 unbiased and perfectly_mappable SNPs: $perfect_n\n\n";
