#!/usr/bin/perl

########################################################################
# 
# 01/16/2012
#
# get_exons_makeCG.pl
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

# this version converts the original exons file to the CG names

use strict;
use warnings;

my$file = $ARGV[0];
my$annot = $ARGV[1];

our($chr,$location,$gene,$range,$strand,$start,$stop,$symbol,$CG);
our(@elements);
our(%names);

open(ANNOT,"$annot") or die "\nError opening $annot\n";
while(<ANNOT>)
{
  chomp;
  if(/^#/){next;}
  else
  {
    @elements = split("\t",$_);
    $CG = $elements[2];
    $symbol = $elements[8];

    unless(!$symbol){$names{$symbol} = $CG;}
  }
}
close ANNOT;

open(FILE,"$file") or die "\nError opening $file\n";
while(<FILE>)
{
  chomp;
  if(/loc=(\S+)\:(\S+)\;\s+name=(\S+)\:+\d+\;/)
  {
    $chr = $1;
    $location = $2;
    $gene = $3;

    if($location =~ /complement\((\d+\.\.\d+)\)/){$range = $1; $strand = "-";} # indicator of strand
    else{$range = $location; $strand = "+";}

    @elements = split(/\.\./,$range);
    $start = $elements[0] - 1; # to make it 0-based for BED format
    $stop = $elements[1];      # already in 1-base

    unless(!$names{$gene}){print "chr$chr\t$start\t$stop\t$names{$gene}\t0\t$strand\n";}
  }
}
close FILE;
