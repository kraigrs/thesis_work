#!/usr/bin/perl

########################################################################
# 
# 01/15/2014
#
# convert_names.pl
#
# Purpose: convert FlyBase FBgn* to CG*
# 
# Input: annotation table, file to be converted 
#
# Output: converted file
#         
# Syntax: perl convert_names <annotations> <list> <output>
#
########################################################################

use strict;
use warnings;

my$annotations = $ARGV[0];
my$list = $ARGV[1];
my$output = $ARGV[2];

my$regulator; my$target; my$score; my$FBgn; my$CG; 
my@elements;
my%FBgn2CG;

open(ANNO,"$annotations") or die "\nError opening $annotations\n";
while(<ANNO>)
{
  chomp;
  unless(/^\#/)
  {
    @elements = split(/\s+/,$_);
    $FBgn = $elements[0];
    $CG = $elements[2];

    $FBgn2CG{$FBgn} = $CG;
  }
}
close ANNO;

open(OUT,"> $output") or die "\nError writing to $output\n";

open(LIST,"$list") or die "\nError opening $list\n";
while(<LIST>)
{
  chomp;
  @elements = split(/\s+/,$_);
  $regulator = $elements[0];
  $target = $elements[1];
  $score = $elements[2];
  
  if($FBgn2CG{$regulator} && $FBgn2CG{$target})
  {
    print OUT "$FBgn2CG{$regulator}\t$FBgn2CG{$target}\t$score\n";
  }
  #else{print "$regulator\n$target\n"}
}
close LIST;

close OUT;
