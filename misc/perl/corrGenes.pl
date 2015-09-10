#!/usr/bin/perl

###################################################################################
#
# 01/12/2010
#
# corrGenes.pl
#
# Purpose: calculate the correlation of log2(S1/S1) between two sets of genes
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

my$spec1 = $ARGV[2];
my$spec2 = $ARGV[3];

my@elements;
my$gene;
my$s1;
my$s2;
my$log2ratio;
my%log2ratios1;
my%spec1_1;
my%spec2_1;

sub log2 
{
  my$n = shift;
  return (log($n)/log(2));
}

open(FILE1,"$file1") or die "\nError opening $file1\n";
while(<FILE1>)
{
  chomp $_;
  @elements = split(/\s/,$_);
  $gene = $elements[0];
  $s1 = $elements[4];
  $s2 = $elements[5];
  
  unless($gene eq "gene" || $gene eq "*")
  {
    unless($s1 == 0 || $s2 == 0)
    {
      $spec1_1{$gene} = $s1;
      $spec2_1{$gene} = $s2;
      $log2ratio = log2($s1/$s2);  
      #print "$log2ratio\n";   
      $log2ratios1{$gene} = $log2ratio;
    }
  }
}
close FILE1;

my%log2ratios2;
my%spec1_2;
my%spec2_2;
  
open(FILE2,"$file2") or die "\nError opening $file2\n";
while(<FILE2>)
{
  chomp $_;
  @elements = split(" ",$_);
  $gene = $elements[0];
  #$s1 = $elements[4];
  #$s2 = $elements[5];
  
  $s1 = $elements[2];
  $s2 = $elements[4];

  #print "$gene\t$s1\t$s2\n";

  unless($gene eq "gene" || $gene eq "*" || $gene eq "Gene")
  {
    unless($s1 == 0 || $s2 == 0)
    {
    $spec1_2{$gene} = $s1;
    $spec2_2{$gene} = $s2;
    $log2ratio = log2($s1/$s2);
    #print "$log2ratio\n";    
    $log2ratios2{$gene} = $log2ratio;
    #print "$log2ratios2{$gene}\t$gene\n";
    }  
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

print "Gene\tAlias\tmel\tsec\tlog2(mel/sec)\tmelHyb\tsecHyb\tlog2(melHyb/secHyb)";

foreach(keys %log2ratios1)
{
  #print "$_\n";
  #print "$log2ratios2{$_}\n";

  if($log2ratios2{$_})
  {
    print "\n$_\t$_\t$spec1_1{$_}\t$spec2_1{$_}\t$log2ratios1{$_}\t$spec1_2{$_}\t$spec2_2{$_}\t$log2ratios2{$_}";
    $n += 1;
    $sumprod += $log2ratios1{$_}*$log2ratios2{$_};
    $sum1 += $log2ratios1{$_};
    $sum2 += $log2ratios2{$_}; 
    $sumsqr1 += $log2ratios1{$_}**2;
    $sumsqr2 += $log2ratios2{$_}**2;
  }
}

my$corr = 0;

$corr = ($n*$sumprod - $sum1*$sum2) / (sqrt($n*$sumsqr1-$sum1**2)*sqrt($n*$sumsqr2-$sum2**2));

print "\nNumber of genes compared for $spec1 X $spec2 hybrid: $n\n";
printf("r = %.3f\n",$corr);
printf("r2 = %.3f\n",$corr**2);

#printf ("\nTime elapsed: %d\n\n",time-$starttime);
