#!/usr/bin/perl

use strict;
use warnings;

my@elements;
my%hash;

my$file = $ARGV[0];
my$out = $ARGV[1];

open(NAMES,"$file") or die "Can't open $file for reading!";
while(<NAMES>)
{
  chomp;
  $hash{$_} = 1;
}
close NAMES;

open(OUT,"> $out") or die "Error writing to $out\n";
foreach(sort keys %hash)
{
  print OUT "$_\n";
}
close OUT;
