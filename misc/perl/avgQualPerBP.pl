#!/usr/bin/perl

########################################################################
# 
# 12/14/2009
#
# avgQualPerBP.pl
#
# Purpose: this script will help in determining a cut-off point for Bowtie
#          by looking at the average read quality per base-pair along
#          a read.
# 
# Input: two files, one for each mate of paired-end reads
#
# Output: two files, one for each mate, containing per base pair position,
#         the average quality over all reads.
#
########################################################################

use strict;
use warnings;

my$start = time;

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

my$mate1 = $ARGV[0];
my$mate2 = $ARGV[1];

my$key1;
my$key2;
my$Ind;
my$numreads1;
my$numreads2;
my$quality;
my%quals1;
my%quals2;
my%qualssqrd1;
my%qualssqrd2;
my@qual;
my$i;
my$avgQual1;
my$avgQual2;
my$stdev1;
my$stdev2;

open(FASTQ1,"$mate1") or die "\nError opening $mate1\n";
$numreads1 = 0;
while(<FASTQ1>)
{
  chomp $_;
  if(/^@(HWI-\S*)\/\d{1}$/)
  {
    $key1 = $1;
    $Ind = 0;
  }
  elsif(/^\+(HWI-\S*)\/\d{1}$/)
  {
    $key2 = $1;
    $Ind = 1;
    $numreads1 += 1;
  }
  else
  {
    if(($Ind == 1) && ($key1 eq $key2))
    {     
      $quality = $_;
      @qual = split("",$_);
      for($i=0;$i<@qual;$i++)
      {
        $quals1{$i+1} += $phred{$qual[$i]};
        $qualssqrd1{$i+1} += $phred{$qual[$i]}**2;
      } 
    }    
  }
}
close FASTQ1;

open(FASTQ2,"$mate2") or die "\nError opening $mate2\n";
$numreads2 = 0;
while(<FASTQ2>)
{
  chomp $_;
  if(/^@(HWI-\S*)\/\d{1}$/)
  {
    $key1 = $1;
    $Ind = 0;
  }
  elsif(/^\+(HWI-\S*)\/\d{1}$/)
  {
    $key2 = $1;
    $Ind = 1;
    $numreads2 += 1;
  }
  else
  {
    if(($Ind == 1) && ($key1 eq $key2))
    {     
      $quality = $_;
      @qual = split("",$_);
      for($i=0;$i<@qual;$i++)
      {
        $quals2{$i+1} += $phred{$qual[$i]};
        $qualssqrd2{$i+1} += $phred{$qual[$i]}**2;
      }  
    }    
  }
}
close FASTQ2;

my$out_dir = $ARGV[2];
my$label = $ARGV[3];

open(OUT1,"> $out_dir/$label.mate1.avgQualPerBP.txt") or die "Error writing to $out_dir/$label.mate1.avgQualPerBP.txt\n";
print OUT1 "base_pair\tavgQual\tstDevQual";
close OUT1;

open(OUT2,"> $out_dir/$label.mate2.avgQualPerBP.txt") or die "Error writing to $out_dir/$label.mate2.avgQualPerBP.txt\n";
print OUT2 "base_pair\tavgQual\tstDevQual";
close OUT2;

foreach(sort {$a <=> $b} keys %quals1)
{
  if($quals2{$_} && $qualssqrd1{$_} && $qualssqrd2{$_} && $numreads1 == $numreads2 )
  {
    $avgQual1 = $quals1{$_} / $numreads1;
    $stdev1 = sqrt((1/($numreads1))*$qualssqrd1{$_} - $avgQual1**2);
    open(OUT1,">> $out_dir/$label.mate1.avgQualPerBP.txt") or die "Error writing to $out_dir/$label.mate1.avgQualPerBP.txt\n";
    print OUT1 "\n$_\t$avgQual1\t$stdev1";
    close OUT1;

    $avgQual2 = $quals2{$_} / $numreads2;
    $stdev2 = sqrt((1/($numreads2))*$qualssqrd2{$_} - $avgQual2**2);
    open(OUT2,">> $out_dir/$label.mate2.avgQualPerBP.txt") or die "Error writing to $out_dir/$label.mate2.avgQualPerBP.txt\n";
    print OUT2 "\n$_\t$avgQual2\t$stdev2";
    close OUT2;
  }
}

printf ("\nTime elapsed: %d\n\n",time-$start);
