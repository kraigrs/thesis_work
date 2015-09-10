#!/usr/bin/perl

######################################################################################
# 
# 01/10/2013
#
# get_5p_3p_regions.pl
#
# Purpose: obtain 5' and 3' extensions of the first and last exons
#
######################################################################################

use strict;
use warnings;

my$ref = $ARGV[0];
my$conversion = $ARGV[1];
my$extend = $ARGV[2];
my$output = $ARGV[3];

my$chr; my$start; my$stop; my$name; my$coords; my$loc; my$exon; my$FBid; my$line;
my$start5p; my$stop5p; my$start3p; my$stop3p;
my@elements; my@meta; my@starts; my@stops;
my%data; my%key;

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
    $exon = $elements[2];
    #print "\n$FBid\t$exon\n\n"; exit;
    unless(!$exon){if($exon =~ /^CG\d+$/){$key{$FBid} = $exon;}}
  }
}
close KEY;

# read in reference sequence
open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  if(/^\>.+loc=(\S+)\:(\S+);\s{1}name=(\S+);\s{1}parent/)
  {
    $chr = $1;
    $coords = $2;
    $name = $3;
    #print "\n$chr\t$coords\t$name\n\n"; exit;

    if($coords =~ /complement\((\d+)\.\.(\d+)\)/)
    {
      $start = $1;
      $stop = $2;
      #print "\n$coords\t$start\t$stop\n\n";
    }
    else
    {
      @meta = split(/\.\./,$coords);
      $start = $meta[0];
      $stop = $meta[1];
      #print "\n$coords\t$start\t$stop\n\n"; exit;
    }
    
    #print "\n$coords\t$start\t$stop\n\n"; exit;

    @meta = split(/\:/,$name);

    if($key{$meta[0]})
    {
      $loc = $key{$meta[0]};
      $exon = $meta[1];

      $data{$loc}{$exon} = [$chr,$start,$stop];
    }
  }
}
close REF;

open(OUT,"> $output") or die "Error writing to $output\n";

foreach $loc (keys %data)
{
  foreach $exon (sort {$a <=> $b} keys %{$data{$loc}})
  {
    $chr = $data{$loc}{$exon}[0];
    push(@starts,$data{$loc}{$exon}[1]);
    push(@stops,$data{$loc}{$exon}[2]);
  }

  $start = shift(@starts);
  $stop = pop(@stops);

  $start5p = $start - $extend - 1;
  $stop5p = $start - 1;

  $start3p = $stop;
  $stop3p = $stop + $extend;

  print OUT "chr$chr\t$start5p\t$stop5p\t$loc\_upstream\n";
  print OUT "chr$chr\t$start3p\t$stop3p\t$loc\_downstream\n";

  @starts = ();
  @stops = ();
}

close OUT;
