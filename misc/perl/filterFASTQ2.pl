#!/usr/bin/perl

########################################################################
# 
# 01/08/2010
#
# filterFASTQ2.pl
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
# Addition: incorporate trimming into the filter
#
########################################################################

use strict;
use warnings;

my$start = time;

my$data_dir = $ARGV[0];
my$fastq = $ARGV[1];

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

my$trim = $ARGV[2];
my$i;

my%seqs1;
my%Ns1;
my%quals1;
my%avgQuals1;

my@seq;
my@qual;

open(FASTQ1,"$data_dir/$fastq\_1_sequence.txt") or die "\nError opening $data_dir/$fastq\_1_sequence.txt\n";
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
      #foreach(@seq){if($_ eq "N"){$Nct += 1;}}
      for($i=0;$i<scalar@seq-$trim;$i++){if($seq[$i] eq "N"){$Nct += 1;}} #trimming
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
        #foreach(@qual){$totqual += $phred{$_};}
        for($i=0;$i<scalar@qual-$trim;$i++){$totqual += $phred{$qual[$i]};} #trimming
        my$avgQual = $totqual / (scalar@qual-$trim);
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

open(FASTQ2,"$data_dir/$fastq\_2_sequence.txt") or die "\nError opening $data_dir/$fastq\_2_sequence.txt\n";
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
      #foreach(@seq){if($_ eq "N"){$Nct += 1;}}
      for($i=0;$i<scalar@seq-$trim;$i++){if($seq[$i] eq "N"){$Nct += 1;}} #trimming
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
        #foreach(@qual){$totqual += $phred{$_};}
        for($i=0;$i<scalar@qual-$trim;$i++){$totqual += $phred{$qual[$i]};} #trimming
        my$avgQual = $totqual / (scalar@qual-$trim);
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
    open(OUT1,">> $data_dir/$fastq\_1_sequence.filtered2.txt") or die "Error writing to $data_dir/$fastq\_1_sequence.filtered2.txt\n";
    print OUT1 "\@$_/1\n$seqs1{$_}\n\+$_/1\n$quals1{$_}\n";
    close OUT1;

    open(OUT2,">> $data_dir/$fastq\_2_sequence.filtered2.txt") or die "Error writing to $data_dir/$fastq\_2_sequence.filtered2.txt\n";
    print OUT2 "\@$_/2\n$seqs2{$_}\n\+$_/2\n$quals2{$_}\n";
    close OUT2;
  }
}

printf ("\nTime elapsed: %d\n\n",time-$start);
