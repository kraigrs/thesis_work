#!/usr/bin/perl

########################################################################
# 
# 01/08/2009
#
# makeFASTQ.pl
#
# Purpose: makes a FASTA file for Bowtie
# 
# Input: a file containing read, tile, lane, mate, qual
#
# Output: FASTQ file
#
########################################################################

use strict;
use warnings;

my$data_dir = $ARGV[0];
my$file = $ARGV[1];
my$suffix = $ARGV[2];

my$name;

open(FILE,"$data_dir/$file") or die "\nError opening $data_dir/$_\n";
while(<FILE>)
{
  chomp $_;
  #print "$_\n";
  if(/^(>HWI\S+)/)
  {
    $name = $1;
  }
  else
  {
    open(OUT,">> $data_dir/$suffix.new.fa") or die "Error writing to $data_dir/$suffix.new.fa\n";
    print OUT "$name/1\n";
    print OUT "$_\n";
    close OUT;
  }
}
close FILE;

