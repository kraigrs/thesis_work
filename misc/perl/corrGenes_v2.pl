#!/usr/bin/perl

###################################################################################
#
# 12/22/2011
#
# corrGenes_v2.pl
#
# Purpose: calculate the correlation of raw counts for allele-specific and total expression
#          with and without accounting for overlapping regions in constitutive exons
# 
# Input: two files, each containing a list of genes and their respective expression
#        counts for each species
#
# Output: the correlation between the log2 ratio of expression
#
###################################################################################

use strict;
use warnings;

#my$starttime = time;

my$file1 = $ARGV[0];
my$file2 = $ARGV[1];
my$comp = $ARGV[2];

our($gene,$s1,$s2,$both);
our(@elements);
our(%file1cts,%file2cts);

open(FILE1,"$file1") or die "\nError opening $file1\n";
while(<FILE1>)
{
  chomp;
  if(/^gene/){next;}
  else
  {
    @elements = split(/\s/,$_);
    $gene = $elements[0];
    $s1 = $elements[1];
    $s2 = $elements[2];
    $both = $elements[3];
  
    $file1cts{$gene} = [$s1,$s2,$both];
  }
}
close FILE1;
  
open(FILE2,"$file2") or die "\nError opening $file2\n";
while(<FILE2>)
{
  chomp;
  if(/^gene/){next;}
  else
  {
    @elements = split(/\s/,$_);
    $gene = $elements[0];
    $s1 = $elements[1];
    $s2 = $elements[2];
    $both = $elements[3];
  
    $file2cts{$gene} = [$s1,$s2,$both];
  }
}
close FILE2;

##################################################
# calculate sample Pearson correlation coefficient

my$n = 0;
my$sumprod = 0;
my$sum1 = 0;
my$sum2 = 0;
my$sumsqr1 = 0;
my$sumsqr2 = 0;
my$corrAS = 0;
my$corrT = 0;

foreach(keys %file1cts)
{
  if($file2cts{$_})
  {
    $n += 2;

    #allele-specific
    $sumprod += $file1cts{$_}[0]*$file2cts{$_}[0];
    $sumprod += $file1cts{$_}[1]*$file2cts{$_}[1];

    $sum1 += $file1cts{$_}[0];
    $sum1 += $file1cts{$_}[1];

    $sum2 += $file2cts{$_}[0]; 
    $sum2 += $file2cts{$_}[1]; 

    $sumsqr1 += $file1cts{$_}[0]**2;
    $sumsqr1 += $file1cts{$_}[1]**2;

    $sumsqr2 += $file2cts{$_}[0]**2;
    $sumsqr2 += $file2cts{$_}[1]**2;
  }
}

print "n=$n\tsumprod=$sumprod\tsum1=$sum1\tsum2=$sum2\tsumsqr1=$sumsqr1\tsumsqr2=$sumsqr2\n"; exit;
$corrAS = ($n*$sumprod - $sum1*$sum2) / (sqrt($n*$sumsqr1-$sum1**2)*sqrt($n*$sumsqr2-$sum2**2));

printf("\nCorrelation of $n allele-specific counts for $comp: r2 = %.3f\n",$corrAS**2);

$n = 0;
$sumprod = 0;
$sum1 = 0;
$sum2 = 0;
$sumsqr1 = 0;
$sumsqr2 = 0;
$corrAS = 0;
$corrT = 0;

foreach(keys %file1cts)
{
  if($file2cts{$_})
  {
    $n += 1;

    #total
    $sumprod += ($file1cts{$_}[0]+$file1cts{$_}[1]+$file1cts{$_}[2])*($file2cts{$_}[0]+$file2cts{$_}[1]+$file2cts{$_}[2]);
    $sum1 += $file1cts{$_}[0]+$file1cts{$_}[1]+$file1cts{$_}[2];
    $sum2 += $file2cts{$_}[0]+$file2cts{$_}[1]+$file2cts{$_}[2]; 
    $sumsqr1 += ($file1cts{$_}[0]+$file1cts{$_}[1]+$file1cts{$_}[2])**2;
    $sumsqr2 += ($file2cts{$_}[0]+$file2cts{$_}[1]+$file2cts{$_}[2])**2;
  }
}

$corrT = ($n*$sumprod - $sum1*$sum2) / (sqrt($n*$sumsqr1-$sum1**2)*sqrt($n*$sumsqr2-$sum2**2));

printf("\nCorrelation of $n total counts for $comp: r2 = %.3f\n",$corrT**2);

sub log2 
{
  my$n = shift;
  return (log($n)/log(2));
}
