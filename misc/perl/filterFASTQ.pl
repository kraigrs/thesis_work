#!/usr/bin/perl

########################################################################
# 
# 10/27/2009
#
# filterFASTQ.pl
#
# Purpose: this filters out mate paired-end reads that each have fewer 
#          than 3 ambiguous bases and avg. quality scores greater or
#          equal to 20
# 
# Input: two files, one for each mate of paired-end reads
#
# Output: two files, one for each mate, each containing the filtered set 
#         of paired-end reads
#
########################################################################

use strict;
use warnings;

my$start = time;

my$data_dir = $ARGV[0];
my$mate1 = $ARGV[1];
my$mate2 = $ARGV[2];

my@list;
my$key;
my%phred;

open(CHARS,"ASCIIchars.txt") or die "\nError opening ASCIIchars.txt\n";
while(<CHARS>)    
{
  chomp $_;
  @list = split(" ",$_);
  $key = $list[0];
  $phred{$key} = $list[2];
  #print "list{$key} --> $phred{$key}\n";
}
close CHARS;

my$key1;
my$key2;
my$Ind;
my$Nct;
my$totqual;

my$sequence;
my$quality;

my%seqs1;
my%Ns1;
my%quals1;
my%avgQuals1;

my@seq;
my@qual;

open(FASTQ1,"$data_dir/$mate1.fq") or die "\nError opening $data_dir/$mate1.fq\n";
while(<FASTQ1>)
{
  chomp $_;
  if(/^@(HWI-\S*)\/\d{1}$/)
  {
    $key1 = $1;
    $Ind = 1;
  }
  elsif(/^\+(HWI-\S*)\/\d{1}$/)
  {
    $key2 = $1;
    $Ind = 0;
  }
  else
  {
    if($Ind == 1)
    {
      $Nct = 0;
      $sequence = $_;
      @seq = split("",$_);
      foreach(@seq){if($_ eq "N"){$Nct += 1;}}
      $seqs1{$key1} = $sequence;
      $Ns1{$key1} = $Nct;
    }
    else
    {
      if($key1 eq $key2)
      {
        $totqual = 0;
        $quality = $_;
        @qual = split("",$_);
        foreach(@qual){$totqual += $phred{$_};}
        my$avgQual = $totqual / scalar @qual;
        $quals1{$key2} = $quality;
        $avgQuals1{$key2} = $avgQual;
      }     
      else{print "\nThe keys do not match! $key1 != $key2\n";}
    }
  }
}
close FASTQ1;

my%seqs2;
my%Ns2;
my%quals2;
my%avgQuals2;

open(FASTQ2,"$data_dir/$mate2.fq") or die "\nError opening $data_dir/$mate2.fq\n";
while(<FASTQ2>)
{
  chomp $_;
  if(/^@(HWI-\S*)\/\d{1}$/)
  {
    $key1 = $1;
    $Ind = 1;
  }
  elsif(/^\+(HWI-\S*)\/\d{1}$/)
  {
    $key2 = $1;
    $Ind = 0;
  }
  else
  {
    if($Ind == 1)
    {
      $Nct = 0;
      $sequence = $_;
      @seq = split("",$_);
      foreach(@seq){if($_ eq "N"){$Nct += 1;}}
      $seqs2{$key1} = $sequence;
      $Ns2{$key1} = $Nct;
    }
    else
    {
      if($key1 eq $key2)
      {
        $totqual = 0;
        $quality = $_;
        @qual = split("",$_);
        foreach(@qual){$totqual += $phred{$_};}
        my$avgQual = $totqual / scalar @qual;
        $quals2{$key2} = $quality;
        $avgQuals2{$key2} = $avgQual;
      }     
      else{print "\nThe keys do not match! $key1 != $key2\n";}
    }
  }
}
close FASTQ2;

foreach(keys %seqs1)
{
  #print "$_\n";
  #print "\@$_\n$seqs1{$_}\n\+$_\n$quals1{$_}\n";
  #print "\@$_\n$seqs2{$_}\n\+$_\n$quals2{$_}\n";
  if($Ns1{$_} < 3 && $avgQuals1{$_} >= 20 && $Ns2{$_} < 3 && $avgQuals2{$_} >= 20)
  {
    open(OUT1,">> $data_dir/$mate1.filtered.fq") or die "Error writing to $data_dir/$mate1.filtered.fq\n";
    print OUT1 "\@$_/1\n$seqs1{$_}\n\+$_/1\n$quals1{$_}\n";
    close OUT1;

    open(OUT2,">> $data_dir/$mate2.filtered.fq") or die "Error writing to $data_dir/$mate2.filtered.fq\n";
    print OUT2 "\@$_/2\n$seqs2{$_}\n\+$_/2\n$quals2{$_}\n";
    close OUT2;
  }
}

printf ("\nTime elapsed: %d\n\n",time-$start);
