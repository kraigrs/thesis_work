#!/usr/bin/perl

use strict;
use warnings;

my$file = $ARGV[0];
my$prefix = $ARGV[1];

my$chr; my$pos; my$refAllele; my$sample; my$i; my$j; my$k; my$columns; my$merge;
my@elements;
my%data;

open(FILE,"$file") or die "Can't open $file for reading!";
while(<FILE>) 
{
  chomp;
  @elements = split(/\s/,$_);
  $chr = $elements[0];
  $pos = $elements[1];
  $refAllele = $elements[2];
  $merge = $pos."_".$refAllele;
  
  $columns = scalar(@elements) - 3;

  for($i = 3; $i < @elements; $i++)
  {
    push @{ $data{$chr}{$merge} }, $elements[$i];
  }
}
close FILE;

for($j = 0; $j < $columns; $j++)
{
  $k = $j + 1;

  open(OUT," > $prefix\_$k.txt") or die "Error writing to $prefix\_$k.txt\n";
  foreach $chr (keys %data)
  {
    foreach $merge (keys %{ $data{$chr} })
    {
      @elements = split(/\_/,$merge);
      $pos = $elements[0];
      $refAllele = $elements[1];
      print OUT "$chr\t$pos\t$refAllele\t$data{$chr}{$merge}[$j]\n";
    }
  }
  close OUT;
}
