#!/usr/bin/perl

########################################################################
# 
# 02/04/2014
#
# get_trim_reads.pl
#
# Purpose: read in two files of unmapped reads to two separate genomes and output a 
#          sequence file of their common reads
# 
# Input: two .un files and the original sequence file
#
# Output: a new .un file representing the set of common reads in each .un file
#
# Syntax: perl get_trim_reads.pl <1.un> <2.un> <original FAST[A/Q] file> <new FAST[A/Q] file> 
#
########################################################################

use strict;
use warnings;

my$file1 = $ARGV[0];
my$file2 = $ARGV[1];
my$seqs = $ARGV[2];
my$output = $ARGV[3];

my$temp; my$name; my$line; my$read; my$p;
my%list;

if($seqs =~ /^(\S+)\.(fastq|fq|fasta|fa)$/){$temp = $1;}
else{$temp = $seqs;}

system "grep '^\@HWI' $file1 | sed 's/^@//' | perl -p -e 's/\\s+\\S+$//\n/' | sort > $file1\.sorted";
system "grep '^\@HWI' $file2 | sed 's/^@//' | perl -p -e 's/\\s+\\S+$//\n/' | sort > $file2\.sorted";
system "comm -12 $file1\.sorted $file2\.sorted > $temp.names";

#system "grep '^\@' $file1 | sed 's/^@//' | perl -p -e 's/\\s+\\S+$//\n/' | sort > $file1\.sorted";
#system "grep '^\@' $file2 | sed 's/^@//' | perl -p -e 's/\\s+\\S+$//\n/' | sort > $file2\.sorted";
#system "comm -12 $file1\.sorted $file2\.sorted > $temp.names";

open(NAMES,"$temp.names") or die "\nError opening $temp.names\n";
while(<NAMES>)
{
  chomp;
  $name = $_;
  $list{$name} = 1;
}
close NAMES;

open(OUT,"> $output") or die "\nError opening $output\n";

$p = 0;
open(SEQS,"$seqs") or die "\nError opening $seqs\n";
while(<SEQS>)
{
  chomp;
  $line = $_;

  if($line =~ /^\@(HWI\S+)/)
  #if($line =~ /^\@(\S+)/)
  {
    $read = $1;
    if($list{$read})
    {
      print OUT "$line\n";
      $p = 1;
    }
    else{$p = 0;}
  }
  elsif($p == 1){print OUT "$line\n";}
}
close SEQS;
close OUT;

system "rm $file1.sorted $file2.sorted $temp.names";
