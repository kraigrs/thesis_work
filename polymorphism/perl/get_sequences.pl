#!/usr/bin/perl

######################################################################################
# 
# 05/02/2012
#
# get_sequences.pl
#
# Purpose: using BED file extract sequences from reference genome
# 
# Input: BED file of exons (or other regions I guess), reference genome (FASTA)        
#
# Output: for each entry in BED file, the actual sequence from the reference
#
# Notes: adapted from the pyroAssay scripts, and greatly improved in speed and memory footprint
#        use this from now on to extract any type of sequences. 
#        could also be used to get single bases, as long as BED file is properly specified
#
######################################################################################

use strict;
use warnings;

my$BED = $ARGV[0];
my$ref = $ARGV[1];
my$output = $ARGV[2];

my$first = 1; my$seq = ""; # without initializing $seq, this is ridiculously slow, 'my' is different from 'our'
our($chr,$start,$stop,$gene,$locus,$begin,$length,$string);
our(@elements,@array);
our(%regs,%genome);

# read in BED file
open(BED,"$BED") or die "\nError opening $BED\n";
while(<BED>)    
{
  chomp;
  @elements = split(/\s+/,$_);
  $chr = $elements[0]; 
  $start = $elements[1];
  $stop = $elements[2];
  $gene = $elements[3];
  $locus = $gene."_".$start."_".$stop;

  #$regs{$chr}{$locus} = [$start,$stop];
  $regs{$chr}{$gene} = [$start,$stop];
}
close BED;

# read in reference sequence
open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  #if(/^\>(\S.chr\S.)$/)
  if(/^\>(\S+)$/)
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
#print "$chr"; exit;

#print "\nDone\n"; exit;

foreach $chr (keys %genome)
{
  $length = length($genome{$chr});
  #print "\nChromosome: $chr\tLength: $length";
}
#print "\n\n";

open(OUT,">$output") or die "Error writing to $output\n";

#print OUT "locus\tsequence";

foreach $chr (keys %regs)
{
  foreach $locus (keys %{$regs{$chr}})
  {
    if($genome{$chr})
    {
      $begin = $regs{$chr}{$locus}[0];
      $length = $regs{$chr}{$locus}[1]-$regs{$chr}{$locus}[0];

      $string = substr $genome{$chr}, $begin, $length;

      print OUT "\n\>$locus\n$string";
    }
  }
}
close OUT;
print "$output\n";
