#!/usr/bin/perl

use strict;
use warnings;

my$fasta = $ARGV[0];
my$length = $ARGV[1];
my$output = $ARGV[2];

my$first = 1; my$seq = "";
my$chr; my$pos; my$exon; my$start; my$stop; my$begin; my$end; my$max;
my%genome;

# read in reference sequence
open(REF,"$fasta") or die "\nError opening $fasta\n";
while(<REF>)    
{
  chomp;
  #if(/^\>(\S.chr\S.)$/)
  if(/^\>(\S+)\|berlin/)
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
$genome{$chr} = $seq;

open(OUT,"> $output") or die "Error writing to $output\n";

foreach $chr (keys %genome)
{
  $begin = 0;
  $max = length($genome{$chr});
    
  while($begin+$length <= $max)
  {
    $end = $begin + $length;
    print OUT "$chr\t$begin\t$end\n";
    $begin += 1;
  }
}
close OUT;


