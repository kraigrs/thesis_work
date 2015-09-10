#!/usr/bin/perl

########################################################################
# 
# 01/13/2013 (will fill in later)
#
########################################################################

use strict;
use warnings;

my$file = $ARGV[0]; #z30

my$file1 = $ARGV[1]; #zhr

my$file2 = $ARGV[2]; #sim
my$geno2 = $ARGV[3];

my$file3 = $ARGV[4]; #sec
my$geno3 = $ARGV[5];

my$out1 = $ARGV[6];
my$out2 = $ARGV[7];

my$first = 1; my$seq = "";
my$chr; my$start; my$stop; my$name; my$base; my$gene; my$locus;
my$major; my$a1; my$a2; my$a3;
my@elements; my@meta;
my%bed; my%bed1; my%bed2; my%bed3; my%genome2; my%genome3;
my%total_locus; my%total_gene; my%specific_locus; my%specific_gene;

#z30
open(BED,"$file") or die "\nError opening $file\n";
while(<BED>)    
{
  chomp;
  @elements = split("\t",$_);
  
  $name = $elements[3];
  $base = $elements[4];

  @meta = split(/\_/,$name);
  if(scalar(@meta) == 5){$locus = $meta[0]."_".$meta[1]."_".$meta[2]."_".$meta[3];}
  else{$locus = $meta[0]."_".$meta[1];}
  $gene = $meta[0];

  $total_locus{$locus} = 0;
  $total_gene{$gene} = 0;

  $specific_locus{$locus} = 0;
  $specific_gene{$gene} = 0;
 
  $bed{$name} = $base;
}
close BED;

#zhr
open(BED,"$file1") or die "\nError opening $file1\n";
while(<BED>)    
{
  chomp;
  @elements = split("\t",$_);
  
  $name = $elements[3];
  $base = $elements[4];
 
  $bed1{$name} = $base;
}
close BED;

#sim
open(BED,"$file2") or die "\nError opening $file2\n";
while(<BED>)    
{
  chomp;
  @elements = split("\t",$_);

  $chr = $elements[0];
  $start = $elements[1];
  $name = $elements[3];
 
  $bed2{$name} = [$chr,$start];
}
close BED;

# read in reference sequence
$first = 1; $seq = "";
open(REF,"$geno2") or die "\nError opening $geno2\n";
while(<REF>)    
{
  chomp;
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

#sec
open(BED,"$file3") or die "\nError opening $file3\n";
while(<BED>)    
{
  chomp;
  @elements = split("\t",$_);

  $chr = $elements[0];
  $start = $elements[1];
  $name = $elements[3];
 
  $bed3{$name} = [$chr,$start];
}
close BED;

# read in reference sequence
$first = 1; $seq = "";
open(REF,"$geno3") or die "\nError opening $geno3\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome3{$chr} = $seq;} # the first time, don't save sequence in hash
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
$genome3{$chr} = $seq; # this needs to be done since last entry will not be stored in hash

# compare alleles

foreach $name (keys %bed)
{
  if($bed1{$name} && $bed2{$name} && $bed3{$name})
  {
    @meta = split(/\_/,$name);
    if(scalar(@meta) == 5){$locus = $meta[0]."_".$meta[1]."_".$meta[2]."_".$meta[3];}
    else{$locus = $meta[0]."_".$meta[1];}
    $gene = $meta[0];

    $total_locus{$locus} += 1;
    $total_gene{$gene} += 1;

    $major = $bed{$name};
    $a1 = $bed1{$name};
    $a2 = substr $genome2{$bed2{$name}[0]},$bed2{$name}[1],1;
    $a3 = substr $genome3{$bed3{$name}[0]},$bed3{$name}[1],1;

    if($major ne $a1 && $a1 eq $a2 && $a2 eq $a3)
    {
      $specific_locus{$locus} += 1;
      $specific_gene{$gene} += 1;
    }
  }
}

open(OUT,">$out1") or die "Error writing to $out1\n";
print OUT "locus\tlineage_specific\ttotal\n";
foreach $locus (keys %total_locus)
{
  print OUT "$locus\t$specific_locus{$locus}\t$total_locus{$locus}\n";
}
close OUT;

open(OUT,">$out2") or die "Error writing to $out2\n";
print OUT "gene\tlineage_specific\ttotal\n";
foreach $gene (keys %total_gene)
{
  print OUT "$gene\t$specific_gene{$gene}\t$total_gene{$gene}\n";
}
close OUT;
