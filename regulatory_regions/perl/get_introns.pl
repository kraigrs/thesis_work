#!/usr/bin/perl

######################################################################################
# 
# 01/10/2013
#
# get_intron_regions.pl
#
# Purpose: get introns from FlyBase db
#
######################################################################################

use strict;
use warnings;

my$ref = $ARGV[0];
my$conversion = $ARGV[1];

my$chr; my$start; my$stop; my$name1; my$name2; my$coords; my$loc; my$intron; my$FBid; my$line; my$left; my$right; my$strand;
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
    $intron = $elements[2];
    unless(!$intron){if($intron =~ /^CG\d+$/){$key{$FBid} = $intron;}}
  }
}
close KEY;

# read in reference sequence
open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  if(/^\>intron\_(\S+)\:(\d+)\_(\S+)\:(\d+).+loc=(\S+)\:(\S+);/)
  {
    $name1 = $1;
    $left = $2;
    $name2 = $3;
    $right = $4;
    $chr = $5;
    $coords = $6;

    #if($name1 ne $name2){print "Check that $name1 and $name2 anre't the same gene!!!\n";}

    if($coords =~ /complement\((\d+)\.\.(\d+)\)/)
    {
      $start = $1 - 1;
      $stop = $2;
      $strand = "-";
      #print "\n$coords\t$start\t$stop\n\n";
    }
    else
    {
      @meta = split(/\.\./,$coords);
      $start = $meta[0] - 1;
      $stop = $meta[1];
      $strand = "+";
      #print "\n$coords\t$start\t$stop\n\n"; exit;
    }
    
    #print "\n$coords\t$start\t$stop\n\n"; exit;

    if($name1 eq $name2)
    {
      if($name1 =~ /CG\d+/)
      { 
        $loc = $name1;
        $intron = "intron_".$loc.":".$left."_".$loc.":".$right;
        #print "\n$intron\n\n"; exit;

        print "chr$chr\t$start\t$stop\t$intron\t0\t$strand\n";
      }
      elsif($key{$name1})
      { 
        $loc = $key{$name1};
        $intron = "intron_".$loc.":".$left."_".$loc.":".$right;
        #print "\n$intron\n\n"; exit;

        print "chr$chr\t$start\t$stop\t$intron\t0\t$strand\n";
      }
    }
  }
}
close REF;
