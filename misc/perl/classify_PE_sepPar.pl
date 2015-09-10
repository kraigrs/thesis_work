#!/usr/bin/perl

###################################################################################
#
# 03/14/2012
#
# classify_PE_sepPar.pl
#
# Purpose: using the merged file from all alignment information, calculate transcript abundances 
# 
# Input: merged file with mate1ref1, mate2ref1
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
#      <merged>   ==> file with mate 1 and mate2 alignments to species1
#      <out_dir>  ==> where to put files
#      <prefix>   ==> what to put at beginning of each file
#
###################################################################################

use strict;
use warnings;

my$const = $ARGV[0];
my$species1 = $ARGV[1]; 
my$merged = $ARGV[2];
my$out_dir = $ARGV[3];
my$prefix = $ARGV[4];

our($line,$gene,$exon,$read,$flag);
our(@elements,@align,@mate1ref1,@mate2ref1,@test);
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
    $gene = $elements[3];
    #$exon = $elements[1] . "," . $elements[2];
    $exon = $elements[1] . "_" . $elements[2];

    $geneHash{$gene}{$species1} = 0; 
    $exonHash{$gene}{$exon}{$species1} = 0;
  }
}
close CONST;
print "\nDone with genes and exons hashes!\n";

# read in the new merged file that has each PE read on a single line and PE-specific alignment info

open(OUT1,">$out_dir\/$prefix\.reads.txt") or die "Error writing to $out_dir\/$prefix\.reads.txt\n";
print OUT1 "read\tmate1_$species1\tmate2_$species1\tcall";

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
  if($elements[2] =~ /No/){@mate2ref1 = ("No","No");}
  else
  {
    @align = split(":",$elements[2]);
    @mate2ref1 = ($align[1],$align[0]);
  }
  
  #print "@mate1ref1\t@mate2ref1\n";
  #print "$mate1ref1[0]\t$mate2ref1[1]\n";

  $flag = &allele_caller(\@mate1ref1,\@mate2ref1);
  print OUT1 "\n$read\t$elements[1]\t$elements[2]\t$flag";
} 
close IN;
close OUT1;

open(OUT2,">$out_dir\/$prefix\.exons.txt") or die "Error writing to $out_dir\/$prefix\.exons.txt\n";
print OUT2 "gene_exon\t$species1";

foreach my$k1 (keys %exonHash)
{
  foreach my$k2 (keys %{$exonHash{$k1}})
  {
    if($exonHash{$k1}{$k2}{$species1} == 0)
    {
      print OUT2 "\n$k1\_$k2\t0";
    }
    else
    {
      print OUT2 "\n$k1\_$k2\t$exonHash{$k1}{$k2}{$species1}";
    }
  }
}
close OUT2;

open(OUT3,">$out_dir\/$prefix\.genes.txt") or die "Error writing to $out_dir\/$prefix\.genes.txt\n";
print OUT3 "gene\t$species1";

foreach(keys %geneHash)
{
  if($geneHash{$_}{$species1} == 0)
  {
    print OUT3 "\n$_\t0";
  }
  else
  {
    print OUT3 "\n$_\t$geneHash{$_}{$species1}";
  }
}
close OUT3;

##### start counting #####
#there should be 4 possibilities for alignments of both mates to both species (2 files, 2 possibilities each, 2^2 = 4)
#these can be taken care of on a case by case basis

sub allele_caller
{
  my($f11,$f21) = @_;
  my$call;
  my@a11 = @{$f11}; my@a21 = @{$f21};

  # mate aligns to both genomes
  if($a11[0] ne "No" && $a21[0] ne "No")
  {
    if($a11[0] eq $a21[0] && $a11[1] eq $a21[1])
    {
      $exonHash{$a11[0]}{$a11[1]}{$species1} += 1; 
      $geneHash{$a11[0]}{$species1} += 1;
      $call = "$species1";
    }
    elsif($a11[0] eq $a21[0])
    {
      $exonHash{$a11[0]}{$a11[1]}{$species1} += 1; 
      $exonHash{$a21[0]}{$a21[1]}{$species1} += 1;
      $geneHash{$a11[0]}{$species1} += 1;
      $call = "$species1";
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
  elsif($a21[0] ne "No")
  {
    $exonHash{$a21[0]}{$a21[1]}{$species1} += 1; 
    $geneHash{$a21[0]}{$species1} += 1;
    $call = "$species1";
  }

  # mate doesnt align
  else{$call = "NA";}

  return $call;
}
