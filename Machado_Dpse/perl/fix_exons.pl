#!/usr/bin/perl

########################################################################
# 
# 11/12/2013
#
# fix_exons.pl
#
# Purpose: take in a VCF file and reference FASTA file and output the alternative reference
# 
# Input: VCF file, original reference, name of alternative reference 
#
# Output: alternative reference in FASTA format
#         
# Syntax: perl fix_exons.pl <FASTA reference> > new.fasta
#
########################################################################

use strict;
use warnings;

my$ref = $ARGV[0];

my$chr; my$coords; my$name; my$start; my$stop; my$exon;
my@meta;

open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  if(/loc=(\S+)\:(\S+);.+parent=(\w+),/)
  {
    $chr = $1;
    $coords = $2;
    $name = $3;

    if($coords =~ /complement\((\d+)\.\.(\d+)\)/)
    {
      $start = $1 - 1;
      $stop = $2;
    }
    else
    {
      @meta = split(/\.\./,$coords);
      $start = $meta[0] - 1;
      $stop = $meta[1];
    }

    $exon = $name."_".$start."_".$stop;
    print "\>$exon\n";
  }
  else{print "$_\n";}
}
close REF;
