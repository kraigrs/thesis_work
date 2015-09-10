#!/usr/bin/perl

########################################################################
# 
# 08/07/2011
#
# split_lane.pl
#
# Purpose: 
# 
# Input: 
#
# Output:
#         
# Syntax: perl split_lane.pl <FASTQ> [lane1, lane2, ...]
#
########################################################################

use strict;
use warnings;

my$file = $ARGV[0];

my$line; my$name; my$seq; my$qual; my$i; my$filename; my$lane; my$string; my$index; my$send; my$outfile;
my@elements; my@meta;
my%lanes; my%outs;

if($file =~ /^(\S+).fastq/){$filename = $1;} 

#for($i=1;$i<@ARGV;$i++)
#{
#  $outfile = "OUT".$ARGV[$i];
#
#  $lanes{$ARGV[$i]} = 1;
#  $outs{$outfile} = 1;
#
#  open(my$outfile,"\> $filename\_$ARGV[$i].fastq") or die "Error writing to $filename\_$ARGV[$i].fastq";
#}

open(OUT2,"\> $filename\_2.fastq") or die "Error writing to $filename\_2.fastq";
open(OUT7,"\> $filename\_7.fastq") or die "Error writing to $filename\_7.fastq";
open(OUT8,"\> $filename\_8.fastq") or die "Error writing to $filename\_8.fastq";

$seq = 0; $qual = 0;
open(IN,"$file") or die "\nError opening $file\n";
while(<IN>)    
{
  chomp;
  $line = $_;

  if(/^\@SRR/)
  {
    $seq = 1;

    @elements = split(/\s+/,$line);
    $name = $elements[1];
    @meta = split(/\:/,$name);
    $lane = $meta[1];

    if($lane == 2){print OUT2 "$line\n";}
    elsif($lane == 7){print OUT7 "$line\n";}
    elsif($lane == 8){print OUT8 "$line\n";}
    else{die "Error with lane numbers!";}
  }
  elsif(/^\+SRR/)
  {
    $qual = 1;

    @elements = split(/\s+/,$line);
    $name = $elements[1];
    @meta = split(/\:/,$name);
    $lane = $meta[1];

    if($lane == 2){print OUT2 "$line\n";}
    elsif($lane == 7){print OUT7 "$line\n";}
    elsif($lane == 8){print OUT8 "$line\n";}
    else{die "Error with lane numbers!";}
  }
  elsif($seq)
  {
    $seq = 0;

    if($lane == 2){print OUT2 "$line\n";}
    elsif($lane == 7){print OUT7 "$line\n";}
    elsif($lane == 8){print OUT8 "$line\n";}
    else{die "Error with lane numbers!";}
  }
  elsif($qual)
  {
    $qual = 0;

    if($lane == 2){print OUT2 "$line\n";}
    elsif($lane == 7){print OUT7 "$line\n";}
    elsif($lane == 8){print OUT8 "$line\n";}
    else{die "Error with lane numbers!";}
  }
}
close IN;

close OUT2; close OUT7; close OUT8;

#foreach(keys %outs){close $_;}
