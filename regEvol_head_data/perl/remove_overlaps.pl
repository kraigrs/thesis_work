#!/usr/bin/perl

use strict;
use warnings;

my$bed = $ARGV[0];
my$out = $ARGV[1]; #list of overlapping regions to remove
#my$outA = $ARGV[1];
#my$outB = $ARGV[2];

our($line,$chrA,$startA,$stopA,$geneA,$chrB,$startB,$stopB,$geneB,$key,$stringA,$stringB);
our(@elements,@A,@B);
our(%overlapA,%overlapB);

open(OUTA,">tempA.bed") or die "\nError writing to tempA.bed\n";
open(OUTB,">tempB.bed") or die "\nError writing to tempB.bed\n";

open(BED,"$bed") or die "\nError opening $bed\n";
while(<BED>)    
{
  chomp;
  $line = $_;

  if(/^track/){next;}
  else
  {
    @elements = split("\t",$line);

    $chrA = $elements[0]; 
    $startA = $elements[1];
    $stopA = $elements[2];
    $geneA = $elements[3];
    $stringA = $chrA."@".$startA."@".$stopA."@".$geneA;

    $chrB = $elements[6]; 
    $startB = $elements[7];
    $stopB = $elements[8];
    $geneB = $elements[9];
    $stringB = $chrB."@".$startB."@".$stopB."@".$geneB;
    
    #print "\n$stringA\t$stringB"; #exit;

    if($overlapB{$stringA} && $overlapA{$stringB})
    {
      #print "\nFound a match";
      @A = split("@",$stringB); #print "\n$overlapA{$stringB}\t$A[0]\n"; exit;
      foreach(@A){print OUTA "$_\t";}
      print OUTA "\n";

      @B = split("@",$stringA);
      foreach(@B){print OUTB "$_\t";}
      print OUTB "\n";
    }
    elsif($stringA ne $stringB)
    {
      $overlapA{$stringA} = 1;
      $overlapB{$stringB} = 1;
    }
  }
}
close BED;
close OUTA; close OUTB;

system "intersectBed -a tempA.bed -b tempB.bed > tempAintempB.bed";
system "intersectBed -a tempB.bed -b tempA.bed > tempBintempA.bed";
system "cat tempAintempB.bed tempBintempA.bed > $out";
