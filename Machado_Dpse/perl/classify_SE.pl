#!/usr/bin/perl

###################################################################################
#
# 02/23/2012
#
# classify_SE.pl
#
# Purpose: using the merged file from all alignment information, calculate transcript abundances 
# 
# Input: merged file with mate1ref1, mate1ref2
#
# Output: a list of genes/exons and their respective match and mismatch counts 
#
# Fixes: this script represents an updated version fixing a memory leak problem. Also, instead of repeatedly
#        building gene and exon hashes, those are made once in the beginning
#
# E.g. perl classify_SE.pl <const> <species1> <species2> <merged> <out_dir> <prefix>
#     
#      <const>    ==> list of constitutive exons
#      <species1> ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <species2> ==> indicate name of 2nd species aligned to
#      <merged>   ==> file with mate 1 alignments to 1st species
#      <out_dir>  ==> where to put files
#      <prefix>   ==> what to put at beginning of each file
#
###################################################################################

use strict;
use warnings;

my$const = $ARGV[0];
my$species1 = $ARGV[1]; 
my$species2 = $ARGV[2];
my$merged = $ARGV[3];
my$out_dir = $ARGV[4];

our($line,$gene,$exon,$read,$flag,$info);
our(@elements,@align,@mate1ref1,@mate1ref2,@test,@meta);
our(%geneHash,%exonHash);

# gene and exon hashes 
open(CONST,"$const") or die "Can't open $const for reading!";
LINE:while(<CONST>) 
{
  chomp;
  if($_ =~ /^track.+/){next LINE;}
  else
  { 
    @elements = split(/\s/,$_);
    #$gene = $elements[3];
    #$exon = $elements[1] . "," . $elements[2];
    #$exon = $elements[1] . "_" . $elements[2];

    $info = $elements[3];
    @meta = split(/\_/,$info);
    $gene = $meta[0];
    $exon = $elements[1]."_".$elements[2];

    $geneHash{$gene}{$species1} = 0; 
    $exonHash{$gene}{$exon}{$species1} = 0;

    $geneHash{$gene}{$species2} = 0; 
    $exonHash{$gene}{$exon}{$species2} = 0;

    $geneHash{$gene}{b} = 0; 
    $exonHash{$gene}{$exon}{b} = 0;
  }
}
close CONST;
print "\nDone with genes and exons hashes!\n";

# read in the new merged file that has each PE read on a single line and PE-specific alignment info

open(OUT1,">$out_dir\.reads.txt") or die "Error writing to $out_dir\.reads.txt\n";
print OUT1 "read\tmate_$species1\tmate_$species2\tcall";

open(IN,"$merged") or die "Can't open $merged for reading!";
while(<IN>) 
{
  chomp;
  @elements = split(/\s+/,$_);
  $read = $elements[0];
  #@test = split("@",$elements[0]); $read = $test[0];

  #print "$read\t$elements[1]\t$elements[2]\t$elements[3]\t$elements[4]\n";

  # mate1ref1
  if($elements[1] =~ /No/){@mate1ref1 = ("No","No");}
  else
  {
    @align = split(":",$elements[1]);
    @mate1ref1 = ($align[1],$align[0]);
  }

  # mate1ref2
  if($elements[2] =~ /No/){@mate1ref2 = ("No","No");}
  else
  {
    @align = split(":",$elements[2]);
    @mate1ref2 = ($align[1],$align[0]);
  }
  
  #print "@mate1ref1\t@mate1ref2\t@mate2ref1\t@mate2ref2\n";
  #print "$mate1ref1[0]\t$mate1ref2[1]\n";

  $flag = &allele_caller(\@mate1ref1,\@mate1ref2);
  print OUT1 "\n$read\t$elements[1]\t$elements[2]\t$flag";
} 
close IN;
close OUT1;

open(OUT2,">$out_dir\.exons.txt") or die "Error writing to $out_dir\.exons.txt\n";
print OUT2 "exon\t$species1\t$species2\tBoth\tAdj$species1\tAdj$species2";

my$madj = 0;
my$sadj = 0;

foreach my$k1 (keys %exonHash)
{
  foreach my$k2 (keys %{$exonHash{$k1}})
  {
    if($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2} == 0)
    {
      print OUT2 "\n$k1\_$k2\t$exonHash{$k1}{$k2}{$species1}\t$exonHash{$k1}{$k2}{$species2}\t$exonHash{$k1}{$k2}{b}\t0\t0";
    }
    else
    {
      $madj = $exonHash{$k1}{$k2}{$species1} + ($exonHash{$k1}{$k2}{$species1} / ($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2})) * $exonHash{$k1}{$k2}{b};  
      $sadj = $exonHash{$k1}{$k2}{$species2} + ($exonHash{$k1}{$k2}{$species2} / ($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2})) * $exonHash{$k1}{$k2}{b};
      print OUT2 "\n$k1\_$k2\t$exonHash{$k1}{$k2}{$species1}\t$exonHash{$k1}{$k2}{$species2}\t$exonHash{$k1}{$k2}{b}\t$madj\t$sadj";
    }
  }
}
close OUT2;

open(OUT3,">$out_dir\.genes.txt") or die "Error writing to $out_dir\.genes.txt\n";
print OUT3 "gene\t$species1\t$species2\tBoth\tAdj_$species1\tAdj_$species2";

$madj = 0;
$sadj = 0;

foreach(keys %geneHash)
{
  if($geneHash{$_}{$species1} + $geneHash{$_}{$species2} == 0)
  {
    print OUT3 "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}\t$geneHash{$_}{b}\t0\t0";
  }
  else
  {
    $madj = $geneHash{$_}{$species1} + ($geneHash{$_}{$species1} / ($geneHash{$_}{$species1} + $geneHash{$_}{$species2})) * $geneHash{$_}{b};
    $sadj = $geneHash{$_}{$species2} + ($geneHash{$_}{$species2} / ($geneHash{$_}{$species1} + $geneHash{$_}{$species2})) * $geneHash{$_}{b};
    print OUT3 "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}\t$geneHash{$_}{b}\t$madj\t$sadj";
  }
}
close OUT3;

##### start counting #####
#there should be 4 possibilities for alignments of both mates to both species (2 files, 2 possibilities each, 2^2 = 4)
#these can be taken care of on a case by case basis

sub allele_caller
{
  my($f11,$f12) = @_;
  my$call;
  my@a11 = @{$f11}; my@a12 = @{$f12};

  # mate aligns to both genomes
  if($a11[0] ne "No" && $a12[0] ne "No")
  {
    if($a11[0] eq $a12[0] && $a11[1] eq $a12[1])
    {
      $exonHash{$a11[0]}{$a11[1]}{b} += 1; 
      $geneHash{$a11[0]}{b} += 1;
      $call = "Both";
    }
    else{$call = "Error";}
  }

  # mate aligns to genome 1 only
  elsif($a11[0] ne "No")
  {
    $exonHash{$a11[0]}{$a11[1]}{$species1} += 1; 
    $geneHash{$a11[0]}{$species1} += 1;
    $call = "$species1";
  }
  
  # mate aligns to genome 2 only
  elsif($a12[0] ne "No")
  {
    $exonHash{$a12[0]}{$a12[1]}{$species2} += 1; 
    $geneHash{$a12[0]}{$species2} += 1;
    $call = "$species2";
  }

  # mate doesnt align
  else{$call = "NA";}

  return $call;
}
