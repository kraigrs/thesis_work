#!/usr/bin/perl

######################################################################################
# 
# 06/04/2012
#
# get_mappability.pl
#
# Purpose: extract mappability measures from GEM output
# 
# Input: GEM mappability output from FASTA genome        
#
# Output: kind of like SAM format, e.g. chr3R 1-based_position mappability
#
######################################################################################

use strict;
use warnings;

my$GEM = $ARGV[0];
my$output = $ARGV[1];

my$encode = 0; my$fasta = 0; my$counter = 0; my$seq = ""; my$first = 1; my$readlen = 0; my$i; my$j;
our($chr,$start,$stop,$gene,$locus,$code,$freq,$begin,$end,$len,$string,$aref,$mappability,$kmer,$SNP);
our(@elements,@array);
our(%map,%genome);

open(GEM,"$GEM") or die "\nError opening $GEM\n";
while(<GEM>)    
{
  chomp;
  #$counter += 1;
  #print "$counter\n";

  if(/READ LENGTH/)
  {
    $readlen = 1;
  }
  elsif($readlen == 1)
  {
    if(/^(\d+)$/)
    {
      $kmer = $1;
      $readlen = 0;
    }
  }

  if(/ENCODING/)
  {
    $encode = 1;
  }
  elsif($encode == 1)
  {
    if(/\'(\S{1})\'\~\[(\d+)\-\d+\]/)
    {
      $code = $1; $freq = $2;
      #print "\nCode: $code\tFreq: $freq\n";
      $map{$code} = $freq;
    }
  }

  if(/^\~(chr\S+)$/)
  {
    unless($first == 1){$genome{$chr} = $seq;}

    $chr = $1;
    $seq = "";
    #print "\nChromosome: $chr\n";

    $fasta = 1;
    $encode = 0;
    $first = 0;
  }
  elsif($fasta == 1)
  {
    $seq = $seq.$_;
    #print "\nSequence: $seq\n";
  }
}
close GEM;
$genome{$chr} = $seq;

#print "\nDone\n"; exit;

#foreach $code (keys %map)
#{
#  print "Code: $code\tFreq: $map{$code}\n";
#}
#exit;

open(OUT,">$output") or die "Error writing to $output\n";
print OUT "chr\tposition\tfrequency";

foreach $chr (keys %genome)
{
  $len = length($genome{$chr});
  print "\nChromosome: $chr\tLength: $len\n";

  for($i=0;$i<$len;$i++)
  {
    $code = substr $genome{$chr}, $i, 1;
    $freq = $map{$code};
    #if(!$map{$code}){print "Chromosome: $chr\tPosition: $i\tCode: $code\n"; exit;}
    if(!$map{$code}){$freq = 0;}
    $j = $i+1;
    print OUT "\n$chr\t$j\t$freq";
  }
}
close OUT;

print "\n";


