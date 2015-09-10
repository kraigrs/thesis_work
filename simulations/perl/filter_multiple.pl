#!/usr/bin/perl

########################################################################
# 
# 02/16/2012
#
# filter_multiple.pl
#
# Purpose: remove multiply mapping reads
#
# Notes: 
# 
# Input: all alignments, multiple alignments
#
# Output: unique alignments  
#
# Syntax: perl confirm_alignments.pl <all> <multiple> <unique>
#
########################################################################

use strict;
use warnings;

my$all = $ARGV[0];
my$multiple = $ARGV[1];
my$unique = $ARGV[2];

our($line,$name);
our(@elements);
our(%m);

open(MULTIPLE,"$multiple") or die "\nError opening $multiple\n";
while(<MULTIPLE>)    
{
  chomp;
  $line = $_;

  if($line =~ /^track/){next;}
  else
  {
    @elements = split("\t",$_);
    $name = $elements[3];
    $m{$name} = 1;
  }
}
close MULTIPLE;

open(UNIQUE,">$unique") or die "\nError writing to $unique\n";

open(ALL,"$all") or die "\nError opening $all\n";
while(<ALL>)    
{
  chomp;
  $line = $_;

  if($line =~ /^track/){next;}
  else
  {
    @elements = split("\t",$_);
    $name = $elements[3];
    unless($m{$name}){print UNIQUE "$line\n";}
  }
}
close ALL;
close UNIQUE;
