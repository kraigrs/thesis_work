#!/usr/bin/perl

########################################################################
# 
# 12/05/2011
#
# filter_expressed.pl
#
# Purpose: filter out exons that were not expressed
#
# Notes: 
# 
# Input: BED file, exons that were actually generated
#
# Output: BED file filtered for expression  
#
# Syntax: perl filter_expressed.pl <BED file> <keep> <output>
#         <BED file>: alignment file in the BED format
#         <keep>:     list of exons to keep
#         <output>:   file to send expressed data to
#
########################################################################

use strict;
use warnings;

my$bed = $ARGV[0];
my$keep = $ARGV[1];
my$output = $ARGV[2];

our($line,$chr,$start,$stop,$gene,$reg);
our(@elements,@meta);
our(%keepHash);

open(KEEP,"$keep") or die "\nError opening $keep\n";
while(<KEEP>)    
{
  chomp;

  @elements = split("\t",$_);
  $gene = $elements[0]; 
  $reg = $elements[1];

  @meta = split(',',$reg);
  $start = $meta[0];
  $stop = $meta[1];
  
  $keepHash{$gene}{$start}{$stop} = 1;
}
close KEEP;

open(OUT,">$output") or die "\nError writing to $output\n";

open(BED,"$bed") or die "\nError opening $bed\n";
while(<BED>)    
{
  chomp;
  if(/^track/){next;}
  else
  {
    $line = $_;

    @elements = split("\t",$_);
    $chr = $elements[0]; 
    $start = $elements[1];
    $stop = $elements[2];
    $gene = $elements[3];

    if($keepHash{$gene}{$start}{$stop})
    {
      print OUT "$line\n";
    }
  }
}
close BED;
close OUT;
