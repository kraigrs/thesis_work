#!/usr/bin/perl

########################################################################
# 
# 01/19/2012
#
# get_transcripts.pl
#
# Purpose: from annotated FlyBase transcripts, define breakpoints, etc. 
# 
# Input: the text or FASTA file from FlyBase 
#
# Output: a BED file containing the exons
#
# Usage:
#
########################################################################

use strict;
use warnings;

my$file = $ARGV[0];

our($chr,$location,$transcript,$range,$comp,$join,$strand,$start,$stop,$count);
our(@elements,@join,@reg);
our(%breaks);

open(FILE,"$file") or die "\nError opening $file\n";
while(<FILE>)
{
  chomp;
  if(/loc=(\S+)\:(\S+)\;.+name=(\S+)\;/)
  {
    $chr = $1;
    $location = $2;
    $transcript = $3;

    #if($location =~ /complement\((\S+)\)/){$range = $1;}
    #elsif($location =~ /join\((\S+)\)/){$range = $1;}
    #else{$range = $location;}

    #print "$chr\t$transcript\t$location\n";

    while($location =~ /(\d+\.\.\d+)+/g)
    {
      $range = $1;
      @reg = split(/\.\./,$range);
      $start = $reg[0] - 1;
      $stop = $reg[1];
      print "chr$chr\t$start\t$stop\t$transcript\n";
    }

    #foreach(@join){print "\t$_";} #print "\n"; #exit;

    #$count = 1;
    #foreach(@join){$breaks{$_} = 1;}
    #foreach(sort {$a <=> $b} keys %breaks)
    #{
    #  $start = $elements[0] - 1; # to make it 0-based for BED format
    #  $stop = $elements[1];      # already in 1-base
    #}

    %breaks = (); @join = (); @elements = ();  @reg = ();
    #print "chr$chr\t$start\t$stop\t$transcript\t0\t$strand\n";
  }
  #elsif(/^\>/){print "$_\n";}
}
close FILE;
#print "Count: $count\n";
