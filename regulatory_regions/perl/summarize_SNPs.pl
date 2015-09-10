#!/usr/bin/perl

########################################################################
# 
# 01/11/2013 (will fill in later)
#
# Example:
#
########################################################################

use strict;
use warnings;

my$SNPs = $ARGV[0];
my$out1 = $ARGV[1];
my$out2 = $ARGV[2];

my$locus; my$gene;
my@elements; my@meta;
my%data1; my%data2;

open(SNP,"$SNPs") or die "\nError opening $SNPs\n";
while(<SNP>)    
{
  chomp;
  @elements = split("\t",$_);
  $locus = $elements[0];
  @meta = split(/\_/,$locus);
  $gene = $meta[0];
  $data1{$locus} += 1;
  $data2{$gene} += 1;
}
close SNP;

open(OUT,"\> $out1") or die "Error opening $out1!";
print OUT "locus\tSNPs\n";
foreach $locus (keys %data1)
{
  print OUT "$locus\t$data1{$locus}\n";
}
close OUT;

open(OUT,"\> $out2") or die "Error opening $out2!";
print OUT "gene\tSNPs\n";
foreach $gene (keys %data2)
{
  print OUT "$gene\t$data2{$gene}\n";
}
close OUT;
