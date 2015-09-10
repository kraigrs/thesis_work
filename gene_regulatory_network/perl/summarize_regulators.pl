#!/usr/bin/perl

########################################################################
# 
# 01/10/2014
#
# summarize_regulators.pl
#
# Purpose: summarize the inferred gene regulatory network from Marbach et al.
# 
# Input: list of regulators, target genes, and score connecting them 
#
# Output: for each regulator --> number of targets and average score for all targets
#         
# Syntax: perl summarize_regulators <list> 
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

  push( @{ $network{$regulator} }, $target);
  push( @{ $support{$regulator} }, $score);
}
close LIST;

$counter = 0;
foreach $regulator (keys %network)
{
  $n = scalar( @{ $network{$regulator} } );
  $sum = 0;

  foreach( @{ $support{$regulator} } ){$sum += $_;}

  $avg = $sum/$n;

  print "$regulator\t$n\t$avg\n";

  $tot += $n;
  $counter += 1;
}

$avg = $tot/$counter;

print "\nAverage number of target genes per regulator: $avg\n\n";
