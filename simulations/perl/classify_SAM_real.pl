#!/usr/bin/perl

###################################################################################
#
# 01/08/2013
#
# classify_SAM_real.pl
#
# Purpose: using existing SNP information, assign reads to an allele
# 
# Input: SAM file that has been aligned to a single reference allowing mismatches 
#
# Output: a list of genes and their respective match and mismatch counts, along with a 
#         summary of all the reads 
#
# E.g. perl classify_SAM.pl 
#
###################################################################################

use strict;
use warnings;

my$ref = $ARGV[0];
my$species1 = $ARGV[1];
my$species2 = $ARGV[2];
my$SNPs = $ARGV[3];
my$sam = $ARGV[4];
my$prefix = $ARGV[5];

my$chr; my$pos; my$base; my$name; my$refct; my$altct; my$start; my$stop; my$flag; my$n; 
my$amb; my$var; my$gene; my$exon; my$gene_exon; my$test1; my$test2; my$test3; my$line;
my$madj; my$sadj;
my@elements; my@read; my@meta;
my%SNPhash; my%refHash; my%geneHash;

my$found = 0; my$seq = "";

#print "\nBuilding reference hash...\n";

open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  #print "Here!\n";
  #if(/^>(chr\S+)$/ && $found == 1)
  if(/^>(\S+)\|berlin/ && $found == 1)
  {
    $refHash{$chr} = $seq;

    $seq = "";
    $chr = $1;
  }
  elsif(/^>(\S+)\|berlin/)
  {
    $seq = "";
    $chr = $1;
    $found = 1;
  }
  elsif($found == 1)
  {
    $seq = $seq.uc($_);
  }
}
$refHash{$chr} = $seq;
close REF;

#foreach(keys %refHash)
#{
#  print "exon: $_\n";
#  $seq = substr($refHash{$_},0,50);
#  print "First 50 bases: $seq\n";
#}
#print "\n";
#exit;

# initialize gene hash 
foreach $exon (keys %refHash)
{
  $geneHash{$exon}{$species1} = 0; 
  $geneHash{$exon}{$species2} = 0; 
  $geneHash{$exon}{b} = 0; 
}

#print "\nBuilding SNP hash...\n";
open(SNPS,"$SNPs") or die "\nError opening $SNPs\n";
while(<SNPS>)
{
  chomp;
  @elements = split("\t",$_);
  $exon = $elements[0];
  $pos = $elements[1];
  $base = $elements[3];

  $SNPhash{$exon}{$pos} = $base;
}
close SNPS;

open(READS,">$prefix\.reads.txt") or die "Error writing to $prefix\.reads.txt\n";
print READS "read\texon\tcall\n";

#print "\nReading SAM file...\n";
open(SAM,"$sam") or die "\nError opening $sam\n";
while (<SAM>) 
{
  chomp;
  $line = $_;

  unless($line =~ /^\@/)
  {
    @elements = split(/\t/,$_);
    unless($elements[2] eq "\*")
    {
      $name = $elements[0];
      $seq = $elements[9];
      $start = $elements[3];
      $stop = $start + length($seq) - 1;
      @read = split("",$seq);
      @meta = split(/\|/,$elements[2]);
      $exon = $meta[0];
      #print "\n$name\t$seq\t$start\t$stop\t$exon\n";
      #exit;

      $refct = 0; $altct = 0;

      foreach $pos (keys %{$SNPhash{$exon}}) 
      {
        if($start <= $pos && $pos <= $stop)
        {
          #print "\n$pos\t$SNPhash{$exon}{$pos}\n\n"; exit;
          #$test1 = $read[$pos-$start];
          #$test1 = $pos-$start;
          #$test2 = substr($refHash{$exon},$pos-1,1);
          #$test3 = $SNPhash{$exon}{$pos};
          #print "\n$name\t$seq\t$start\t$stop\t$exon\n";
          #print "$pos\t$test1\t$test2\t$test3\n";
          #exit;

          if($read[$pos-$start] eq substr($refHash{$exon},$pos-1,1)){$refct += 1;}
          elsif($read[$pos-$start] eq $SNPhash{$exon}{$pos}){$altct += 1;}
	}
      }
    
      if($refct > 0 && $altct == 0)
      {
        $flag = $species1;
        $geneHash{$exon}{$species1} += 1; 
      }
      elsif($refct == 0 && $altct > 0)
      {
        $flag = $species2;
        $geneHash{$exon}{$species2} += 1; 
      }
      elsif($refct == 0 && $altct == 0)
      {
        $flag = "Both";
        $geneHash{$exon}{b} += 1; 
      }
      else{$flag = "error";}
      #print "$flag\n";
      print READS "$name\t$exon\t$flag\n";
    }
  }
}
close READS;

open(EXONS,">$prefix\.exons.txt") or die "Error writing to $prefix\.exons.txt\n";
print EXONS "exon\t$species1\t$species2\tBoth\tAdj_$species1\tAdj_$species2";

$madj = 0;
$sadj = 0;

foreach(keys %geneHash)
{
  if($geneHash{$_}{$species1} + $geneHash{$_}{$species2} == 0)
  {
    print EXONS "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}\t$geneHash{$_}{b}\t0\t0";
  }
  else
  {
    $madj = $geneHash{$_}{$species1} + ($geneHash{$_}{$species1} / ($geneHash{$_}{$species1} + $geneHash{$_}{$species2})) * $geneHash{$_}{b};
    $sadj = $geneHash{$_}{$species2} + ($geneHash{$_}{$species2} / ($geneHash{$_}{$species1} + $geneHash{$_}{$species2})) * $geneHash{$_}{b};
    print EXONS "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}\t$geneHash{$_}{b}\t$madj\t$sadj";
  }
}
close EXONS;
