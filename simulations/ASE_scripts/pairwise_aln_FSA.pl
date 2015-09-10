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

my$ref1 = $ARGV[0];
my$ref2 = $ARGV[1];

my$loc; my$seq; my$seq1; my$seq2; my$first; my$length1; my$length2; my$name1; my$name2; my$line;
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

$length1 = 0;
foreach(keys %genome1)
{
  $length1 += length($genome1{$_});
}
print "\n$ref1: $length1\n";

$length2 = 0;
foreach(keys %genome2)
{
  $length2 += length($genome2{$_});
}
print "\n$ref2: $length2\n\n";
exit;

#my$exon = "F51776_SI";
#print "\>$exon:berlin\n$genome1{$exon}\n";
#print "\>$exon:c1674\n$genome2{$exon}\n";
#exit;

if($ref1 =~ /^(\S+).fasta/){$out1 = $1."_fsa.fasta";}
#print "$out1\n";
if($ref2 =~ /^(\S+).fasta/){$out2 = $1."_fsa.fasta";}
#print "$out2\n";

open(OUT1,">$out1") or die "Error opening $out1!";
open(OUT2,">$out2") or die "Error opening $out2!";

# compare common features between genomes

foreach $loc (keys %features)
{
  if($genome1{$loc} && $genome2{$loc})
  {
    system "echo \"\>$loc\_berlin\n$genome1{$loc}\" \> target";
    system "echo \"\>$loc\_c1674\n$genome2{$loc}\" \> query";

    system "fsa target query --stockholm > alignment";

    open(ALN,"alignment") or die "\nError opening alignment\n";
    while(<ALN>)    
    {
      chomp;
      $line = $_;

      if(/_berlin\s+/)
      {
        @elements = split(/\s+/,$line);
        $name1 = $elements[0];
        $seq1 = $elements[1];
        #print "\>$name\n$seq\n";
        $length1 = length($seq1);
      }
      elsif(/_c1674\s+/)
      {
        @elements = split(/\s+/,$line);
        $name2 = $elements[0];
        $seq2 = $elements[1];
        #print "\>$name\n$seq\n";
        $length2 = length($seq2);
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

system "rm target query alignment";
