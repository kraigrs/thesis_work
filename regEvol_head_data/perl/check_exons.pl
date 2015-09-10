#!/usr/bin/perl

use strict;
use warnings;

my$file = $ARGV[0];

my$line; my$chr; my$start; my$stop; my$gene; my$loc; my$strand; my$exon; my$key; my$ref; my$count;
my@elements; my@stuff; my@meta; my@genes;
my%check;

open(BED,"$file") or die "\nError opening $file\n";
while(<BED>)    
{
  chomp;
  $line = $_;

  if(/^track/){next;}
  else
  {
    @elements = split("\t",$line);

    $chr = $elements[0]; 
    $start = $elements[1];
    $stop = $elements[2];
    $loc = $elements[3];
    $strand = $elements[4];

    @stuff = split(/\;/,$loc);

    if(scalar(@stuff) < 2){next;}
    else
    {
      foreach(@stuff)
      {
        @meta = split(/\_/,$_);
        push @genes, $meta[0];
        $ref = $meta[0];
      }
      foreach(@genes)
      {
        #print "$loc\n$ref\n$_\n"; exit;
        if($_ ne $ref)
        {
          print "$loc\n";
          $check{$_} = 1;
          $check{$ref} = 1;
        }
      }
      @genes = ();
    }
  }
}
close BED;

$count = keys %check;
print "\n$count\n\n";
