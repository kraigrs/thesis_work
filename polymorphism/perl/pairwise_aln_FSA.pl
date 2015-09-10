#!/usr/bin/perl

########################################################################
# 
# 07/30/2011
#
# pairwise_aln_FSA.pl
#
# Purpose: compare common fasta features between two files
# 
# Input: two fasta files with some common features (chromosomes, exons, genes, etc.) 
#
# Output: lengths of common features, # that differ, etc.
#         
# Syntax: perl compare_fasta.pl <FASTA1> <FASTA2>
#
########################################################################

use strict;
use warnings;

my$spec1 = $ARGV[0];
my$spec2 = $ARGV[1];
my$ref1 = $ARGV[2];
my$ref2 = $ARGV[3];

my$loc; my$seq; my$seq1; my$seq2; my$first; my$length1; my$length2; my$name1; my$name2; my$line; my$found;
my$out1; my$out2;
my@elements; my@seq1; my@seq2;
my%genome1; my%genome2; my%features;

# read in reference sequences

$first = 1; $seq = "";
open(REF,"$ref1") or die "\nError opening $ref1\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome1{$loc} = $seq;} # the first time, don't save sequence in hash
    #@elements = split(/\|/,$1);
    $loc = $1;
    $features{$loc} = 1;
    $seq = "";
    $first = 0;
    #print "\nChromosome: $loc\n";
  }
  else
  {
    $seq = $seq.uc($_);
    #print "\nSequence: $seq\n";
  }
}
close REF;
$genome1{$loc} = $seq; # this needs to be done since last entry will not be stored in hash
#print "\n\nReference1 stored\n";

$first = 1; $seq = "";
open(REF,"$ref2") or die "\nError opening $ref2\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome2{$loc} = $seq;} # the first time, don't save sequence in hash
    #@elements = split(/\|/,$1);
    $loc = $1;
    $features{$loc} = 1;
    $seq = "";
    $first = 0;
    #print "\nChromosome: $loc\n";
  }
  else
  {
    $seq = $seq.uc($_);
    #print "\nSequence: $seq\n";
  }
}
close REF;
$genome2{$loc} = $seq; # this needs to be done since last entry will not be stored in hash
#print "\n\nReference2 stored\n";

if($ref1 =~ /^(\S+)\.[fasta|fa]/){$out1 = $1."\.$spec1\_$spec2\_fsa.fa";}
#print "$out1\n";
if($ref2 =~ /^(\S+)\.[fasta|fa]/){$out2 = $1."\.$spec1\_$spec2\_fsa.fa";}
#print "$out2\n";

open(OUT1,">$out1") or die "Error opening $out1!";
open(OUT2,">$out2") or die "Error opening $out2!";

# compare common features between genomes

foreach $loc (keys %features)
{
  if($genome1{$loc} && $genome2{$loc})
  {
    #system "echo \"\>$loc\_berlin\n$genome1{$loc}\" \> target";
    #system "echo \"\>$loc\_c1674\n$genome2{$loc}\" \> query";

    system "echo \"\>$loc\_$spec1\n$genome1{$loc}\" \> target_$spec1\_$spec2";
    system "echo \"\>$loc\_$spec2\n$genome2{$loc}\" \> query_$spec1\_$spec2";

    system "fsa target_$spec1\_$spec2 query_$spec1\_$spec2 --stockholm > alignment_$spec1\_$spec2";
    #exit;

    $found = 0;
    open(ALN,"alignment_$spec1\_$spec2") or die "\nError opening alignment_$spec1\_$spec2\n";
    while(<ALN>)    
    {
      chomp;
      $line = $_;

      if($found == 0)
      {
        unless(/^[\#|\/\/]/)
        {
          @elements = split(/\s+/,$line);
          $name1 = $elements[0];
          $seq1 = $elements[1];
          #print "\>$name1\n$seq1\n"; exit;
          #print "\>$name1\n"; exit;
          $length1 = length($seq1);
          #print "length: $length1\n"; exit;
          $found = 1;
        }
      }
      elsif($found == 1)
      {
        unless(/^[\#|\/\/]/)
	{
          @elements = split(/\s+/,$line);
          $name2 = $elements[0];
          $seq2 = $elements[1];
          #print "\>$name\n$seq\n";
          $length2 = length($seq2);
	}
      }
    }
    close ALN;
    #print "$name1\t$name2\n";

    #$length1 = length($genome1{$loc});
    #$length2 = length($genome2{$loc});

    #print "$loc\t$length1\t$length2\n";

    #if($length1 != $length2)
    #{
    #  system "echo \"\>$loc\_berlin\n$genome1{$loc}\" \> target";
    #  system "echo \"\>$loc\_c1674\n$genome2{$loc}\" \> query";
    #  exit;
    #}

    if(!$length1 || !$length2)
    {
      print "$loc\n";
    }
    else
    {
      print OUT1 "\>$loc\n$seq1\n";
      print OUT2 "\>$loc\n$seq2\n";
    }
  }
}

#system "rm target query alignment";
