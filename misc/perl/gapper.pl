#!/usr/bin/perl

########################################################################
# 
# 01/28/2009
#
# gapper.pl
#
# Purpose: reformat some text
# 
# Input: some sort of sequence information
#
# Output: file in right format
#
########################################################################

use strict;
use warnings;

my$file = $ARGV[0];
#my$output = $ARGV[1];

my@list;

my$obs;
my$length;
my$numline;
my$taxon; 
my$read;
my%seqs;
my@seq;

open(FILE,"$file") or die "\nError opening $file\n";
while(<FILE>)
{
  chomp $_;
  #print "$_\n";
  if(/^(\d+)\s+(\d+)$/)
  {
    $obs = $1; $length = $2;
  }
  elsif(/^(\s+.+)/)
  {
    $numline = $1;
  }
  else
  {
    if(/^(\w+)\s+(.+)$/)
    {
      $taxon = $1; $read = $2;
      $seqs{$taxon} = $read;   
    }
  }
}
close FILE;

#print "$obs\t$length";
#print "\n$numline";

my$i;
my@nums;
my@regs;

foreach(keys %seqs)
{
  #print "\n$_\t$seqs{$_}";
  @seq = split(/\s+/,$seqs{$_});
  for($i=0;$i<@seq;$i++)
  {
    if($seq[$i] eq "-"){push @nums,$i}
  }
}
@nums = sort {$a<=>$b} @nums;
#foreach(@nums){print "$_\n";}

for($i=0;$i<@nums;$i++)
{
  if($nums[$i+1] != $nums[$i]){push @regs,$nums[$i];}
}
@regs = sort {$a<=>$b} @regs;
foreach(@regs){print "$_\n";}
    
#open(OUT,">> $output") or die "Error writing to $output\n";
#print OUT "";
#close OUT;
