#!/usr/bin/perl

########################################################################
# 
# 06/14/2011
#
# read_trimmer.pl
#
# Purpose: take a FASTA or FASTQ file and trims the reads (for programs that can't do this themselves... MOSAIK!!!) 
# 
# Input: FASTA or FASTQ file 
#
# Output: the same FASTA or FASTQ file with reads trimmed (either from the 5' or 3' end, or both)
#
# Usage: read_trimmer.pl <fasta or fastq file> <# reads to trim from 5' end> <# reads to trim from 3' end> <output>
#
########################################################################

use strict;
use warnings;

my$file = $ARGV[0];
my$trim5 = $ARGV[1];
my$trim3 = $ARGV[2];
my$out = $ARGV[3];

my$length;
my$trimmed;

open(OUT,"> $out") or die "Error writing to $out\n";

if($file =~ /^\S+\.(fastq|fq)$/)
{
  open(IN,"$file") or die "\nError opening $file\n";
  while(<IN>)
  {
    chomp;
    if(/^\@.+/){print OUT "$_\n";}
    elsif(/^\+.+/){print OUT "$_\n";}
    else
    {
      $length = length($_);

      if($length-$trim5-$trim3 <= 0){die "\nWARNING: the reads should be greater than 1 bp after trimming!\n\nKilling read_trimmer.pl\n";}

      $trimmed = substr($_,$trim5,$length-$trim5-$trim3); # accepts 5' only, 3' only, and both trimming
      print OUT "$trimmed\n";
    }
  }
}
elsif($file =~ /^\S+\.(fasta|fa)$/)
{
  open(IN,"$file") or die "\nError opening $file\n";
  while(<IN>)
  {
    chomp;
    if(/^\>.+/){print OUT "$_\n";}
    else
    {
      $length = length($_);

      if($length-$trim5-$trim3 <= 0){die "\nWARNING: the reads should be greater than 1 bp after trimming!\n\nKilling read_trimmer.pl\n\n";}

      $trimmed = substr($_,$trim5,$length-$trim5-$trim3); # accepts 5' only, 3' only, and both trimming
      print OUT "$trimmed\n";
    }
  }
}
close IN;
close OUT;
