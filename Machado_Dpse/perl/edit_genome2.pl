#!/usr/bin/perl

########################################################################
# 
# 10/04/2011
#
# edit_genome.pl
#
# Purpose: take in a reference genome and return the edited alternate genome and ambiguous genome
# 
# Input: SNP file, original reference, directory for new references 
#
# Output: outputs bad SNPs (SNP showing inconsistency with reference)
#         
# Syntax: perl edit_genome2.pl <genome> <genotypes> <name of alternate genome> > <BED file of SNPs> 
#
########################################################################

use strict;
use warnings;
use Switch;

my$ref = $ARGV[0];
my$line_genos = $ARGV[1];
my$alt_genome = $ARGV[2];

my$chr; my$pos; my$seq=""; my$i; my$base; my$j=0; my$k=0; my$first=1; my$len; my$ref_allele; my$one;
my@elements;
my%genome; my%genos;

# read in reference sequence
open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  if(/^\>(chr\S+)$/)
  {
    unless($first == 1){$genome{$chr} = $seq;} # the first time, don't save sequence in hash
    $chr = $1;
    $seq = "";
    $first = 0;
    #print "\nChromosome: $chr\n";
  }
  else
  {
    $seq = $seq.uc($_);
    #print "\nSequence: $seq\n";
  }
}
close REF;
$genome{$chr} = $seq; # this needs to be done since last entry will not be stored in hash

#print "\n\nReference stored\n";

open(GENO,"$line_genos") or die "\nError opening $line_genos\n";
while(<GENO>)
{
  chomp;  
  @elements = split(/\s+/,$_);
  $chr = $elements[0];
  $pos = $elements[1];
  $base = $elements[3];
  $genos{$chr}{$pos} = $base;   
}
close GENO;

#specify alternate genomes
open(ALT,">$alt_genome") or die "\nError opening $alt_genome\n";

foreach $chr (keys %genome)
{
  print ALT "\>$chr\n";
  $len = length($genome{$chr});

  for($i=0;$i<$len;$i++) 
  {
    $ref_allele = substr $genome{$chr}, $i, 1;

    if($ref_allele eq "N")
    {
      print ALT "N";
      $j += 1;
    }
    elsif($genos{$chr}{$i})
    {
      if($ref_allele ne $genos{$chr}{$i})
      {
        $k += 1;
        $one = $i+1; 
        print "$chr\t$i\t$one\t$ref_allele\t$genos{$chr}{$i}\n";
        print ALT "$genos{$chr}{$i}";
      }
      else
      {
        print ALT "$ref_allele";
      }
      $j += 1;
    }
    else
    {
      print ALT "$ref_allele";
      $j += 1;  
    }
  
    if($j == 50){print ALT "\n"; $j=0;}
  }

  unless($j == 0){print ALT "\n";}
  $j = 0;
}
close ALT;

#print "\nNumber mismatches between reference and DGRP SNP calls: $k\n\n";
