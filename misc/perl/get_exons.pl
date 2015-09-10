#!/usr/bin/perl

########################################################################
# 
# 01/16/2012
#
# get_exons.pl
#
# Purpose: from annotated FlyBase exons, create a BED file of these exons 
# 
# Input: the text or FASTA file from FlyBase 
#
# Output: a BED file containing the exons
#
# Usage:
#
########################################################################

use strict;
use warnings;

my$file = $ARGV[0];

our($chr,$location,$gene,$range,$strand,$start,$stop);
our(@elements);

my$count = 0;

open(FILE,"$file") or die "\nError opening $file\n";
while(<FILE>)
{
  chomp;
  if(/loc=(\S+)\:(\S+)\;\s+name=(\S+\:+\d+)\;/)
  {
    #$count++;
    $chr = $1;
    $location = $2;
    $gene = $3;

    if($location =~ /complement\((\d+\.\.\d+)\)/){$range = $1; $strand = "-";} # indicator of strand
    else{$range = $location; $strand = "+";}

    @elements = split(/\.\./,$range);
    $start = $elements[0] - 1; # to make it 0-based for BED format
    $stop = $elements[1];      # already in 1-base

    print "chr$chr\t$start\t$stop\t$gene\t0\t$strand\n";
  }
  #elsif(/^\>/){print "$_\n";}
}
close FILE;
#print "Count: $count\n";
