#!/usr/bin/perl

###################################################################################
#
# 12/20/2013
#
# reformatBed.pl
#
# Purpose: make BED file amenable to characterize and classify scripts     
#
###################################################################################

use strict;
use warnings;

my$bed = $ARGV[0];

my$gene; my$read; my$start; my$stop;
my@elements; my@meta;

open(BED,"$bed") or die "Can't open $bed for reading!";
while(<BED>) 
{
  chomp;
  unless(/^#/)
  { 
    @elements = split(/\s/,$_);
    @meta = split(/\_/,$elements[0]);
    $read = $elements[3];
    $gene = $meta[0];
    $start = $meta[1];
    $stop = $meta[2];

    print "$gene\t$start\t$stop\t$read\n";
  }
}
close BED;
