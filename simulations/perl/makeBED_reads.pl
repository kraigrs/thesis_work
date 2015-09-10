#!/usr/bin/perl

###################################################################################
#
# 02/13/2011
#
# makeBED_reads.pl
#
# Purpose: make BED files for reads that failed to align to genome
# 
# Input: unaligned reads (FASTA or FASTQ) 
#
# Output: same file but in BED format, using the read annotations 
#
# Usage: perl makeBED_reads.pl <reads> <BED file>
#        <reads>    ==> file containing FASTA or FASTQ reads 
#        <BED file> ==> name of BED file to send reads to
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
  my$reads = $ARGV[0];
  my$BEDfile = $ARGV[1];

  our($name,$chr,$start,$stop);
  our(@elements);

  open(BED,">$BEDfile") or die "Can't open $BEDfile for reading!";

  open(READS,"$reads") or die "Can't open $reads for reading!";
  while(<READS>) 
  {
    chomp;
    if($_ =~ /^\>(HWI\S+)/)
    {
      $name = $1;
      @elements = split(/_/,$name);
      $chr = $elements[1];
      $start = $elements[3];
      $stop = $elements[4];
      print BED "chr$chr\t$start\t$stop\t$name\n";

      #print "$name\n@elements\n";
      #exit;
    }
  }
  close READS;
  close BED;
}
