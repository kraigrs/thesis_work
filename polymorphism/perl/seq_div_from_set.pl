#!/usr/bin/perl

###################################################################################
#
# 02/22/2013
#
# seq_div_from_set.pl
#
###################################################################################

use strict;
use warnings;

my$list = $ARGV[0];
my$SNPs = $ARGV[1];
my$lengths = $ARGV[2];

my$gene; my$exon; my$start; my$stop; my$length; my$name; my$total; my$site;
my@elements; my@meta;
my%list; my%coding; my%geneHash;

open(LIST,"$list") or die "Can't open $list for reading!";
while(<LIST>) 
{
  chomp;
  $list{$_} = 1;
}
close LIST;

open(LENGTH,"$lengths") or die "Can't open $lengths for reading!";
while(<LENGTH>) 
{
  chomp;

  @elements = split(/\s/,$_);
  $gene = $elements[0];
  $length = $elements[0];

  $coding{$gene} = $length;
}
close LENGTH;

open(SNP,"$SNPs") or die "Can't open $SNPs for reading!";
while(<SNP>) 
{
  chomp;

  @elements = split(/\s/,$_);
  $name = $elements[0];
  @meta = split(/\_/,$name);
  $gene = $meta[0];

  $geneHash{$gene} += 1;
}
close SNP;

$total = 0; $sites = 0;
foreach $gene (keys %list)
{
  if($coding{$gene} && $geneHash{$gene})
  {
    $total += $coding{$gene};
    $sites += $geneHash{$gene};
  }
}

print "Total coding bases: $total\n";
print "Total differentiating sites: $sites\n";
