#!/usr/bin/perl

######################################################################################
# 
# 02/12/2013
#
# get_intergenic_regions.pl
#
# Purpose: get intergenic regions from FlyBase annotations
#
######################################################################################

use strict;
use warnings;

my$ref = $ARGV[0];
my$conversion = $ARGV[1];

my$chr; my$start; my$stop; my$name; my$coords; my$locus; my$inter; my$FBid; my$line; my$left; my$right;
my@elements; my@meta;
my%key;

# read in key
open(KEY,"$conversion") or die "\nError opening $conversion\n";
while(<KEY>)    
{
  chomp;
  $line = $_;
  unless($line =~ /^\#/)
  {
    #print "$line\n";
    @elements = split("\t",$line);
    $FBid = $elements[0];
    $locus = $elements[2];
    $key{$FBid} = $locus;

    #unless(!$locus){if($locus =~ /^CG\d+$/){$key{$FBid} = $locus;}}
  }
}
close KEY;

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
      if($key{$left} && $key{$right} && ($key{$left} =~ /CG/ || $key{$right} =~ /CG/) )
      {
        print "chr$chr\t$start\t$stop\t$key{$left}\_intergenic\_$key{$right}\t0\t\+\n";
      }
    }
    elsif($right)
    {
      if($key{$right} && $key{$right} =~ /CG/)
      {
        print "chr$chr\t$start\t$stop\tintergenic\_$key{$right}\t0\t\+\n";
      }
    }
    elsif($left)
    {
      if($key{$left} && $key{$left} =~ /CG/)
      {
        print "chr$chr\t$start\t$stop\t$key{$left}\_intergenic\t0\t\+\n";
      }
    }
  }
}
close REF;
