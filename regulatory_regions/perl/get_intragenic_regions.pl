#!/usr/bin/perl

######################################################################################
# 
# 02/20/2013
#
# get_intragenic_regions.pl
#
# Purpose: intergenic regions for genes that are contained within introns of other genes
#
######################################################################################

use strict;
use warnings;

my$ref = $ARGV[0];
my$conversion = $ARGV[1];

my$chr; my$start; my$stop; my$name; my$coords; my$locus; my$intron; my$line; my$left; my$right; my$length;
my@elements; my@meta;
my%key;

# read in key
open(KEY,"$conversion") or die "\nError opening $conversion\n";
while(<KEY>)    
{
  chomp;
  $line = $_;

  @elements = split(/\s+/,$line);
  $locus = $elements[3];
  $intron = $elements[9];

  $key{$intron}{$locus} = 1;
}
close KEY;

foreach $intron (keys %key)
{
  $length = scalar(keys %{ $key{$intron} });
  #if($length ==17){foreach $locus (keys %{ $key{$intron} }){print "$locus\n";};}
  print "$length\n";
}
exit;

# read in reference sequence
open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S*)\_intergenic\_(\S+)\s{1}type.+loc=(\S+)\:(\S+);\s{1}/)
  {
    $left = $1;
    $right = $2;
    $chr = $3;
    $coords = $4;

    if($coords =~ /complement\((\d+)\.\.(\d+)\)/)
    {
      $start = $1 - 1;
      $stop = $2;
      #print "\n$coords\t$start\t$stop\n\n";
    }
    else
    {
      @meta = split(/\.\./,$coords);
      $start = $meta[0] - 1;
      $stop = $meta[1];
      #print "\n$coords\t$start\t$stop\n\n"; exit;
    }
    
    #print "\n$coords\t$start\t$stop\n\n"; exit;

    if($left && $right)
    {
      if($key{$left} && $key{$right})
      {
        print "$chr\t$start\t$stop\t$key{$left}\_intergenic\_$key{$right}\t0\t\+\n";
      }
    }
    elsif($right)
    {
      if($key{$right})
      {
        print "$chr\t$start\t$stop\tintergenic\_$key{$right}\t0\t\+\n";
      }
    }
    elsif($left)
    {
      if($key{$left})
      {
        print "$chr\t$start\t$stop\t$key{$left}\_intergenic\t0\t\+\n";
      }
    }
  }
}
close REF;
