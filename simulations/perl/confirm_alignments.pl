#!/usr/bin/perl

########################################################################
# 
# 12/05/2011
#
# confirm_alignments.pl
#
# Purpose: for reads generated from genome, make sure that read metadata matches target alignment
#
# Notes: 
# 
# Input: BED file
#
# Output: summary statistics: # reads that don't match, output reads that don't match to another file  
#
# Syntax: perl confirm_alignments.pl <BED file> <output>
#         <BED file>: alignment file in the BED format
#         <output>:   a list of reads that fail to match their metadata to target alignment
#
########################################################################

use strict;
use warnings;

my$bed = $ARGV[0];
#my$output = $ARGV[1];

our($line,$chr,$start,$stop,$name,$false,$total,$metaChr,$metaStart,$metaStop,$label,$n);
our(@elements,@meta);
our(%dups);

if($bed =~ /^(.*)\.bed$/){$label = $1."\.false.bed";}
$false = 0; $total = 0;

#open(OUT,">$output") or die "\nError writing to $output\n";
open(OUT,">$label") or die "\nError writing to $label\n";

$n = 0;

open(BED,"$bed") or die "\nError opening $bed\n";
while(<BED>)    
{
  chomp;
  if(/^track/){next;}
  else
  {
    $total += 1;
    $line = $_;

    @elements = split("\t",$_);
    $chr = $elements[0]; 
    $start = $elements[1];
    $stop = $elements[2];
    $name = $elements[3];

    if($dups{$name}){$n += 1;}
    else{$dups{$name} = 1;}

    @meta = split('_',$name);
    $metaChr = 'chr'.$meta[1];
    $metaStart = $meta[3];
    $metaStop = $meta[4];

    #print "$chr = $metaChr\n$start = $metaStart\n$stop = $metaStop\n"; exit;

    if($chr ne $metaChr || $start != $metaStart || $stop != $metaStop)
    {
      $false += 1;
      #print "$chr = $metaChr\n$start = $metaStart\n$stop = $metaStop\n"; exit;
      print OUT "$line\n";
    }

  }
}
close BED;
close OUT;

print "\nAlignments that appear to be false: $false\nTotal successful alignments: $total\nNumber of non-unique alignments: $n\n\n";
