#!/usr/bin/perl

use strict;
use warnings;

my$fasta = $ARGV[0];

my$chr; my$first; my$seq; my$length; my$total; my$total_chr; my$GC; my$GC_chr;
my$i; my$base; my$percent;
my@elements;
my%genome;

# initialize variables

$first = 1;
$seq = "";

# read in fasta sequence(s)

open(FASTA,"$fasta") or die "\nError opening $fasta\n";
while(<FASTA>)    
{
  chomp;
  if(/^\>(\S+)\s*/)
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
close FASTA;
$genome{$chr} = $seq; # this needs to be done since last entry will not be stored in hash

#print "Done reading in FASTA"; 

print "feature\tpercent_GC";

$total = 0; $GC = 0;
foreach $chr (keys %genome)
{
  $length = length($genome{$chr});
  $total_chr = 0;
  $GC_chr = 0;

  for($i = 0; $i < $length; $i++)
  {
    $base = substr($genome{$chr},$i,1);
    if($base =~ /[ACTGactg]/)
    {
      $total += 1;
      $total_chr += 1;
      if($base =~ /[GCgc]/)
      {
        $GC += 1;
        $GC_chr += 1;
      }
    }
  }
  $percent = $GC_chr/$total_chr;
  print "\n$chr\t$percent";
}
$percent = $GC/$total;
print "\nTotal\t$percent\n";
