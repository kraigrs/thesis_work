#!/usr/bin/perl

########################################################################
# 
# 01/3/2011
#
# convert2FASTA.pl
#
# Purpose: converts FASTQ (or mixture of FASTQ and FASTA) to FASTA
# 
# Input: a file containing read, tile, lane, mate, qual, etc.
#
# Output: FASTA file
#
########################################################################

use strict;
use warnings;

my$file = $ARGV[0];

my$name;
my$found = 0;
my$prefix;

if($file =~ /^(.+)\.f\w+/){$prefix = $1;}
print "\n$prefix\n\n";

open(OUT,">>$prefix.fa") or die "Error writing to $prefix.fa\n";

open(FILE,"$file") or die "\nError opening $file\n";
while(<FILE>)
{
  chomp $_;
  #print "$_\n";

  if($found)
  {
    $found = 0;
    print OUT "\>$name\n";
    print OUT "$_\n";
  }
  else
  {
    if(/^\>(HWI\S+)/)
    {
      $found = 1;
      $name = $1;
    }
    elsif(/^\@(HWI\S+)/)
    {
      $found = 1;
      $name = $1;
    }
  }
}
close OUT;
close FILE;

