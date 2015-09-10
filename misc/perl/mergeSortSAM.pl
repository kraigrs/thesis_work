#!/usr/bin/perl

###################################################################################
#
# 02/08/2010
#
# mergeSortSAM.pl
#
# Purpose: merge and sort 2 .sam files (ideally from two mates from the same experiment)
# 
# Input: the two .sam files (mate 1 and mate 2). The alignments to a particular gene
#        will be checked to make sure they are the same. If only 1 mate aligned, then
#        that will be counted as a single read. Otherwise, both mate paired end reads
#        must align to the same gene and will be counted as a single read
#
# Output: a concatenated .sam file containing all reads from the mate pairs (if applicable) 
#
# Example: perl mergeSortSAM.pl ../mel_sec_data/Hyb_mate1.dmel-gene.sam ../mel_sec_data/Hyb_mate2.dmel-gene.sam ../mel_sec_data/Hyb.dmel-gene.sam ../mel_sec_data/Hyb_pairs.dmel-gene.sam ../mel_sec_data/Hyb_singles.dmel-gene.sam > ../mel_sec_data/mergeSortSAM.txt &
#
###################################################################################

use strict;
use warnings;

my$SAM_file1 = $ARGV[0];
my$SAM_file2 = $ARGV[1];
my$SAMout = $ARGV[2];
my$matepairs = $ARGV[3];
my$singlereads = $ARGV[4]; 

my$mate1found = 0;
my$mate2found = 0;
my$read1;
my$read2;
my$gene1;
my$gene2;

my@elements;

my$sam1;
my$sam2;

#print "\n\nParsing first mates file...\n";
open(SAM1,"$SAM_file1") or die "\nError opening $SAM_file1\n";
while(<SAM1>)
{
  chomp $_;  
  # make the header
  if(/^@/)
  {
    open(OUT1,">> $SAMout") or die "Error writing to $SAMout\n";
    print OUT1 "$_\n";
    close OUT1;
  }
  elsif(/^HWI.+\/1\s+/)
  {
    #print "\nMate 1 found!\n";
    $mate1found = 1;
    @elements = split(/\s/,$_);
    #print "$elements[0]";
    #print "$_\n";
    if($elements[0] =~ /^(HWI.+)\/\d{1}/){$read1 = $1;}
    #print "mate 1 found: $read1\n";
    $gene1 = $elements[2]; 
    $sam1 = $_;

    open(SAM2,"$SAM_file2") or die "\nError opening $SAM_file2\n";
    PARSE: while(<SAM2>)
    {
      #print "Parsing second mates file";
      chomp $_;   
      if(/^HWI.+\/2\s+/)
      {
        @elements = split(/\s/,$_);
        #print "$elements[0]";
        #print "$_\n";
        if($elements[0] =~ /^(HWI.+)\/\d{1}/){$read2 = $1;}
        #print "mate 2 found: $read2\n";
        $gene2 = $elements[2];
        $sam2 = $_;
        if(($read1 eq $read2) && ($gene1 eq $gene2) && ($gene1 ne "*"))
	{
          open(OUT1,">> $SAMout") or die "Error writing to $SAMout\n";
          print OUT1 "$sam1\n$sam2\n";
          close OUT1;

          open(OUT2,">> $matepairs") or die "Error writing to $matepairs\n";
          print OUT2 "$sam1\n$sam2\n";
          close OUT2; 

          #print "Match!\n$sam1\n$sam2\n";
          last PARSE;
	}
        elsif(($read1 eq $read2) && ($gene1 ne "*"))
	{
          open(OUT1,">> $SAMout") or die "Error writing to $SAMout\n";
          print OUT1 "$sam1\n";
          close OUT1; 

          open(OUT3,">> $singlereads") or die "Error writing to $singlereads\n";
          print OUT3 "$sam1\n";
          close OUT3;

          #print "Mate2 didn't align!\n$sam1\n$sam2\n";
          last PARSE;
	}
        elsif(($read1 eq $read2) && ($gene2 ne "*"))
	{
          open(OUT1,">> $SAMout") or die "Error writing to $SAMout\n";
          print OUT1 "$sam2\n";
          close OUT1; 

          open(OUT3,">> $singlereads") or die "Error writing to $singlereads\n";
          print OUT3 "$sam2\n";
          close OUT3; 

          #print "Mate1 didn't align!\n$sam1\n$sam2\n"; 
          last PARSE;
	}
        elsif(($read1 eq $read2) && ($gene1 eq $gene2) && ($gene1 eq "*"))
	{
          #print "Neither mate aligned!\n$sam1\n$sam2\n";
          last PARSE;
        }
      }
    }
    close SAM2;
  }
  #print "Closed mate 2 file...\n\n";
}
close SAM1;
#print "Finished!\n"
