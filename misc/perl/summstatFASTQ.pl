#!/usr/bin/perl

########################################################################
# 
# 10/09/2009
#
# summstatFASTQ.pl
#
# Purpose: collect summary statistics for FASTQ reads, including number
#          of mismatches and average quality score
# 
# Input: a file with FASTQ reads
#
# Output: a table with, for each read, the number of mismatches and 
#         the quality of the read (average of the base qualities?)
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

#opendir(DIR,$data_dir) || die "\nCan't open directory $data_dir for reading\n";
#my@files = grep{/^s/} grep{/.txt$/} readdir(DIR);
#close(DIR);

my$key1;
my$key2;
my$Ind;
my$Nct;
my$totqual;

my@seq;
my@qual;

#foreach(@files)
#{
  #open(OUT,"> $_\.summary.txt") or die "Error writing to $_\.summary\n";
  open(OUT,"> $fastq.summary.txt") or die "Error writing to $fastq.summary.txt\n";
  print OUT "read\tNprop\tavgQual";
  close OUT;

  #open(FASTQ,"$data_dir/$_") or die "\nError opening $data_dir/$_\n";
  open(FASTQ,"$data_dir/$fastq.txt") or die "\nError opening $data_dir/$fastq.txt\n";
  while(<FASTQ>)
  {
    chomp $_;
    #print "$_\n";  
    if(/^@(HWI-\S*)/)
    {
      $key1 = $1;
      $Ind = 1;
    }
    elsif(/^\+(HWI-\S*)/)
    {
      $key2 = $1;
      $Ind = 0;
    }
    else
    {
      if($Ind == 1)
      {
        $Nct = 0;
        @seq = split("",$_);
        foreach(@seq)
        {
          if($_ eq "N"){$Nct += 1;}
        }
        #open(OUT,"> $_\.summary.txt") or die "Error writing to $_\.summary\n";
        open(OUT,">> $fastq.summary.txt") or die "Error writing to $fastq.summary.txt\n";
        print OUT "\n$key1\t$Nct";
        close OUT;
      }
      else
      {
        if($key1 eq $key2)
	{
          $totqual = 0;
          @qual = split("",$_);
          foreach(@qual)
          {
            $totqual += $phred{$_};
          }
          my$avgQual = $totqual / scalar @qual;
          #open(OUT,"> $_\.summary.txt") or die "Error writing to $_\.summary\n";
          open(OUT,">> $fastq.summary.txt") or die "Error writing to $fastq.summary.txt\n";
          print OUT "\t$avgQual";
          close OUT;  
        }
        else{print "Keys do not match!";}
      }
    }
  }
  close OUT;
  close FASTQ;
#}

#system "say \"Ready to analyze!\"";

printf ("\nTime elapsed: %d\n\n",time-$start);
