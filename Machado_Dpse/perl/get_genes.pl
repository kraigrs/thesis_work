#!/usr/bin/perl

######################################################################################
# 
# 02/15/2013
#
# get_genes.pl
#
# Purpose: get genes from FlyBase annotation
#
# perl get_genes.pl <annotation set> <conversion>
#
######################################################################################

use strict;
use warnings;

my$ref = $ARGV[0];
#my$conversion = $ARGV[1];

my$chr; my$start; my$stop; my$name; my$coords; my$loc; my$gene; my$FBid; my$line; my$strand;
my@elements; my@meta;
my%key;

# read in key
#open(KEY,"$conversion") or die "\nError opening $conversion\n";
#while(<KEY>)    
#{
#  chomp;
#  $line = $_;
#  unless($line =~ /^\#/)
#  {
#    #print "$line\n";
#    @elements = split("\t",$line);
#    $FBid = $elements[0];
#    $gene = $elements[2];
#    unless(!$gene){if($gene =~ /^CG\d+$/){$key{$FBid} = $gene;}}
#  }
#}
#close KEY;

# read in reference sequence
open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)\s+.+loc=(\S+)\:(\S+);/)
  {

    $chr = $2;
    $coords = $3;
    #$name = $1;
    $loc = $1;

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

    #if($key{$name})
    #{
    #  $loc = $key{$name};
    #  #print "\n$intron\n\n"; exit;

      print "$chr\t$start\t$stop\t$loc\t0\t$strand\n";
    #}
  }
}
close REF;
