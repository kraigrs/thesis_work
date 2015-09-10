#!/usr/bin/perl

########################################################################
# 
# 01/15/2014
#
# summarize_targets.pl
#
# Purpose: summarize the inferred gene regulatory network from Marbach et al.
# 
# Input: list of regulators, target genes, and score connecting them 
#
# Output: for each gene --> number of regulators and average score for all connections
#         
# Syntax: perl summarize_targets <list> 
#
########################################################################

use strict;
use warnings;

my$list = $ARGV[0];

my$regulator; my$target; my$score; my$n; my$sum; my$avg; my$tot; my$counter;
my@elements;
my%network; my%support;

open(LIST,"$list") or die "\nError opening $list\n";
while(<LIST>)
{
  chomp;
  @elements = split(/\s+/,$_);
  $regulator = $elements[0];
  $target = $elements[1];
  $score = $elements[2];
  #print "$regulator\t$target\t$score\n";

  push( @{ $network{$target} }, $regulator);
  push( @{ $support{$target} }, $score);
}
close LIST;

$counter = 0;
foreach $target (keys %network)
{
  $n = scalar( @{ $network{$target} } );
  $sum = 0;

  foreach( @{ $support{$target} } ){$sum += $_;}

  $avg = $sum/$n;

  print "$target\t$n\t$avg\n";

  $tot += $n;
  $counter += 1;
}

$avg = $tot/$counter;

print "\nAverage number of regulators for each target gene: $avg\n\n";
