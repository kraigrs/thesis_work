#!/usr/bin/perl

######################################################################################
# 
# 05/02/2012
#
# find_indels.pl
#
# Purpose: determine whether or not an indel neighbors a differentiating site
#
######################################################################################

use strict;
use warnings;

my$ref1 = $ARGV[0];
my$ref2 = $ARGV[1];
my$SNPs = $ARGV[2];
my$length = $ARGV[3];
my$out = $ARGV[4];

my$first = 1; my$seq = "";
my$chr; my$pos; my$seq1; my$seq2; my$indel1; my$indel2;
my@elements;
my%genome1; my%genome2; my%SNPhash;


# read in reference sequence
$first = 1; $seq = ""; $chr = "";
open(REF,"$ref1") or die "\nError opening $ref1\n";
while(<REF>)    
{
  chomp;
  #if(/^\>(\S.chr\S.)$/)
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome1{$chr} = $seq;} # the first time, don't save sequence in hash
    $chr = $1;
    $seq = "";
    $first = 0;
    #print "\nChromosome: $chr\n";
  }
  else
  {
    $seq = $seq.uc($_);
    #print "\nSequence: $seq\n";
  }
}
close REF;
$genome1{$chr} = $seq; # this needs to be done since last entry will not be stored in hash


# read in alternative sequence
$first = 1; $seq = ""; $chr = "";
open(REF,"$ref2") or die "\nError opening $ref2\n";
while(<REF>)    
{
  chomp;
  #if(/^\>(\S.chr\S.)$/)
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome2{$chr} = $seq;} # the first time, don't save sequence in hash
    $chr = $1;
    $seq = "";
    $first = 0;
    #print "\nChromosome: $chr\n";
  }
  else
  {
    $seq = $seq.uc($_);
    #print "\nSequence: $seq\n";
  }
}
close REF;
$genome2{$chr} = $seq; # this needs to be done since last entry will not be stored in hash

open(SNP,"$SNPs") or die "\nError opening $SNPs\n";
while(<SNP>)
{
  chomp;
  @elements = split("\t",$_);

  $chr = $elements[0];
  $pos = $elements[1]; # 1-based position in SNP BED file since 1-based in pileup

  $SNPhash{$chr}{$pos} = 1;
}
close SNP;

open(OUT,"\> $out") or die "Error opening $out!";
print OUT "chr\tpos\tref_indel\talt_indel\n";

foreach $chr (keys %SNPhash)
{
  foreach $pos (keys %{ $SNPhash{$chr} })
  {
    $seq1 = substr $genome1{$chr}, $pos-$length, 2*$length-1;
    $seq2 = substr $genome2{$chr}, $pos-$length, 2*$length-1;

    #print "$seq1\n$seq2\n"; exit;

    if($seq1 =~ /\-/){$indel1 = 1;}
    else{$indel1 = 0;}

    if($seq2 =~ /\-/){$indel2 = 1;}
    else{$indel2 = 0;}

    print OUT "$chr\t$pos\t$indel1\t$indel2\n";
  }
}
close OUT;
