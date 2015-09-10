#!/usr/bin/perl

########################################################################
# 
# 01/08/2009
#
# makeFASTQ.pl
#
# Purpose: makes a FASTQ file for Bowtie
# 
# Input: a file containing read, tile, lane, mate, qual
#
# Output: FASTQ file
#
########################################################################

use strict;
use warnings;

my$start = time;

my$data_dir = $ARGV[0];

opendir(DIR,$data_dir) || die "\nCan't open directory $data_dir for reading\n";
my@files = grep{/^s/} grep{/.txt$/} readdir(DIR);
#my@files = grep{/^s_\d{1}_\d{1}_test/} readdir(DIR);
close(DIR);

my$suffix;
my@list;
my$tech;
my$lane;
my$tile;
my$num1;
my$num2;
my$mate;
my$read;
my$qual;

#foreach(@files){print "$_\n";}
foreach(@files)
{
  open(FILE,"$data_dir/$_") or die "\nError opening $data_dir/$_\n";
  if($_ =~ /^(s_\d{1}_\d{1})/){$suffix = $1;}
  #print "$suffix\n";
  while(<FILE>)
  {
    chomp $_;
    #print "$_\n";
    @list = split(" ",$_); 
    $tech = $list[0];
    $lane = $list[1];
    $tile = $list[2];
    $num1 = $list[3];
    $num2 = $list[4];
    $mate = $list[5];
    $read = $list[6];
    $qual = $list[7];

    #print ">$tech\_lane\_$lane\_tile\_$tile\_$num1\_$num2";
    #print "\n$read";
    #print "\n+$tech\_lane\_$lane\_tile\_$tile\_$num1\_$num2";
    #print "\n$qual\n";
    
    open(OUT,">> $data_dir/$suffix\_sequence.txt") or die "Error writing to $data_dir/$suffix\_sequence.txt\n";
    print OUT ">$tech\_lane\_$lane\_tile\_$tile\_mate\_$mate\_$num1\_$num2";
    print OUT "\n$read";
    print OUT "\n+$tech\_lane\_$lane\_tile\_$tile\_mate\_$mate\_$num1\_$num2";
    print OUT "\n$qual\n";
    close OUT;
  }
  close FILE;
}
