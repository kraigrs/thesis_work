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

opendir(DIR,$data_dir) || die "\nCan't open directory $data_dir for reading\n";
my@files = grep{/^Hyb/} grep{/.fq$/} readdir(DIR);
close(DIR);

my$suffix;
my$name;
my$ind = 0;

#foreach(@files){print "$_\n";}
foreach(@files)
{
  open(FILE,"$data_dir/$_") or die "\nError opening $data_dir/$_\n";
  if($_ =~ /^(Hyb\S+).fq$/){$suffix = $1;}
  #print "$suffix\n";
  while(<FILE>)
  {
    chomp $_;
    #print "$_\n";
    if(/^@(HWI\S+)/)
    {
      $name = $1;
      $ind = 1;
    }
    elsif(/^\+HWI\S+/)
    {
      $ind = 0;
    }
    else
    {
      if($ind == 1) 
      {
        open(OUT,">> $data_dir/$suffix.fa") or die "Error writing to $data_dir/$suffix.fa\n";
        print OUT ">$name\n";
        print OUT "$_\n";
        close OUT;
      }
    }
  }
  close FILE;
}
