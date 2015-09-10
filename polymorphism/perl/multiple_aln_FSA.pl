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
# Syntax: perl multiple_aln_FSA.pl <FASTA1> <FASTA2> [<FASTA3> ...]
#
########################################################################

use strict;
use warnings;

my$label1 = $ARGV[0];
my$label2 = $ARGV[1];
my$label3 = $ARGV[2];
my$label4 = $ARGV[3];
my$label5 = $ARGV[4];

my$ref1 = $ARGV[5];
my$ref2 = $ARGV[6];
my$ref3 = $ARGV[7];
my$ref4 = $ARGV[8];
my$ref5 = $ARGV[9];

my$loc; my$seq; my$first; my$line;
my$seq1; my$seq2; my$seq3; my$seq4; my$seq5;
my$name1; my$name2; my$name3; my$name4; my$name5;
my$length1; my$length2; my$length3; my$length4; my$length5;
my$out1; my$out2; my$out3; my$out4; my$out5;
my@elements; my@seq1; my@seq2;
my%genome1; my%genome2; my%genome3; my%genome4; my%genome5; my%features;

# read in reference sequences

$first = 1; $seq = "";
open(REF,"$ref1") or die "\nError opening $ref1\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome1{$loc} = $seq;} # the first time, don't save sequence in hash
    @elements = split(/\|/,$1);
    $loc = $elements[0];
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
    @elements = split(/\|/,$1);
    $loc = $elements[0];
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

$first = 1; $seq = "";
open(REF,"$ref3") or die "\nError opening $ref3\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome3{$loc} = $seq;} # the first time, don't save sequence in hash
    @elements = split(/\|/,$1);
    $loc = $elements[0];
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
$genome3{$loc} = $seq; # this needs to be done since last entry will not be stored in hash
#print "\n\nReference2 stored\n";

$first = 1; $seq = "";
open(REF,"$ref4") or die "\nError opening $ref4\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome4{$loc} = $seq;} # the first time, don't save sequence in hash
    @elements = split(/\|/,$1);
    $loc = $elements[0];
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
$genome4{$loc} = $seq; # this needs to be done since last entry will not be stored in hash
#print "\n\nReference2 stored\n";

$first = 1; $seq = "";
open(REF,"$ref5") or die "\nError opening $ref5\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome5{$loc} = $seq;} # the first time, don't save sequence in hash
    @elements = split(/\|/,$1);
    $loc = $elements[0];
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
$genome5{$loc} = $seq; # this needs to be done since last entry will not be stored in hash
#print "\n\nReference2 stored\n";

if($ref1 =~ /^(\S+).fa/){$out1 = $1."\.$label1\_$label2\_$label3\_$label4\_$label5\_fsa.fa";}
#print "$out1\n";
if($ref2 =~ /^(\S+).fa/){$out2 = $1."\.$label1\_$label2\_$label3\_$label4\_$label5\_fsa.fa";}
#print "$out2\n";
if($ref3 =~ /^(\S+).fa/){$out3 = $1."\.$label1\_$label2\_$label3\_$label4\_$label5\_fsa.fa";}
#print "$out3\n";
if($ref4 =~ /^(\S+).fa/){$out4 = $1."\.$label1\_$label2\_$label3\_$label4\_$label5\_fsa.fa";}
#print "$out4\n";
if($ref5 =~ /^(\S+).fa/){$out5 = $1."\.$label1\_$label2\_$label3\_$label4\_$label5\_fsa.fa";}
#print "$out5\n";

open(OUT1,">$out1") or die "Error opening $out1!";
open(OUT2,">$out2") or die "Error opening $out2!";
open(OUT3,">$out3") or die "Error opening $out3!";
open(OUT4,">$out4") or die "Error opening $out4!";
open(OUT5,">$out5") or die "Error opening $out5!";

# compare common features between genomes

foreach $loc (keys %features)
{
  if($genome1{$loc} && $genome2{$loc} && $genome3{$loc}&& $genome4{$loc}&& $genome5{$loc})
  {
    system "echo \"\>$loc\_ref1\n$genome1{$loc}\n".
           "\>$loc\_ref2\n$genome2{$loc}\n".
           "\>$loc\_ref3\n$genome3{$loc}\n".
           "\>$loc\_ref4\n$genome4{$loc}\n".
           "\>$loc\_ref5\n$genome5{$loc}\n".
           "\" \> seqs";

    system "fsa --stockholm seqs > alignment";

    open(ALN,"alignment") or die "\nError opening alignment\n";
    while(<ALN>)    
    {
      chomp;
      $line = $_;

      if(/_ref1\s+/)
      {
        @elements = split(/\s+/,$line);
        $name1 = $elements[0];
        $seq1 = $elements[1];
        #print "\>$name\n$seq\n";
        $length1 = length($seq1);
      }
      elsif(/_ref2\s+/)
      {
        @elements = split(/\s+/,$line);
        $name2 = $elements[0];
        $seq2 = $elements[1];
        #print "\>$name\n$seq\n";
        $length2 = length($seq2);
      }
      elsif(/_ref3\s+/)
      {
        @elements = split(/\s+/,$line);
        $name3 = $elements[0];
        $seq3 = $elements[1];
        #print "\>$name\n$seq\n";
        $length3 = length($seq3);
      }
      elsif(/_ref4\s+/)
      {
        @elements = split(/\s+/,$line);
        $name4 = $elements[0];
        $seq4 = $elements[1];
        #print "\>$name\n$seq\n";
        $length4 = length($seq4);
      }
      elsif(/_ref5\s+/)
      {
        @elements = split(/\s+/,$line);
        $name5 = $elements[0];
        $seq5 = $elements[1];
        #print "\>$name\n$seq\n";
        $length5 = length($seq5);
      }
    }
    close ALN;

    #$length1 = length($genome1{$loc});
    #$length2 = length($genome2{$loc});

    #print "$loc\t$length1\t$length2\n";

    #if($length1 != $length2)
    #{
    #  system "echo \"\>$loc\_berlin\n$genome1{$loc}\" \> target";
    #  system "echo \"\>$loc\_c1674\n$genome2{$loc}\" \> query";
    #  exit;
    #}

    if(!$length1 || !$length2 || !$length3 || !$length4 || !$length5)
    {
      #print "$loc\n";
      next;
    }
    else
    {
      print OUT1 "\>$loc\n$seq1\n";
      print OUT2 "\>$loc\n$seq2\n";
      print OUT3 "\>$loc\n$seq3\n";
      print OUT4 "\>$loc\n$seq4\n";
      print OUT5 "\>$loc\n$seq5\n";
    }
  }
}

system "rm seqs alignment";

print "$out1\n$out2\n$out3\n$out4\n$out5\n";
