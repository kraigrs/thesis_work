#!/usr/bin/perl

######################################################################################
# 
# 01/17/2012
#
# merge_SNPs_mapp.pl
#
######################################################################################

use strict;
use warnings;

my$SNPs = $ARGV[0];
my$mapp = $ARGV[1];
my$read_length = $ARGV[2];
my$mis = $ARGV[3];
my$output = $ARGV[4];

my$chr; my$pos; my$base_ref; my$base_alt; my$mappability;
my$i; my$j; my$k; my$right; my$left; my$steps_right; my$steps_left;
my@elements; my@positions;
my%SNPhash; my%mappHash; my%neighbors;

open(SNPS,"$SNPs") or die "\nError opening $SNPs\n";
while(<SNPS>)
{
  chomp;
  @elements = split("\t",$_);
  $chr = $elements[0];
  $pos = $elements[1];
  $base_ref = $elements[2];
  $base_alt = $elements[3];

  $SNPhash{$chr}{$pos} = [$base_ref,$base_alt];
}
close SNPS;

open(MAP,"$mapp") or die "\nError opening $mapp\n";
while(<MAP>)
{
  chomp;
  unless(/^locus/)
  {
    @elements = split("\t",$_);
    $chr = $elements[0];
    $pos = $elements[1];
    $mappability = $elements[2]/$elements[3];

    $mappHash{$chr}{$pos} = $mappability;
  }
}
close MAP;

# algorithm to compute the number of neighbors
foreach $chr (keys %SNPhash)
{
  foreach $pos (sort {$a <=> $b} keys %{$SNPhash{$chr}})
  {
    push @positions, $pos;
  }

  if(scalar(@positions) == 1){$neighbors{$chr}{$positions[0]} = 0;}
  else
  {
    for($i = 0; $i < @positions; $i++)
    {
      #print "$i\n";
      if($i == 0)
      {
        $steps_right = 0;
        $k = $i + 1;
        $right = $positions[$k] - $positions[$i];

        while($right < $read_length)
        {
          $steps_right += 1;
          if($k == scalar(@positions) - 1){last;}
          $k += 1;
          $right = $positions[$k] - $positions[$i];  
        }

        $neighbors{$chr}{$positions[$i]} = $steps_right;
      }
      elsif($i == scalar(@positions) - 1)
      {
        $steps_left = 0;
        $j = $i - 1;
        $left = $positions[$i] - $positions[$j];

        while($left < $read_length)
        {
          $steps_left += 1;
          if($j == 0){last;}
          $j -= 1;
          $left = $positions[$i] - $positions[$j];  
        }

        $neighbors{$chr}{$positions[$i]} = $steps_left;       
      }
      else
      {
        $steps_right = 0;
        $steps_left = 0;

        $k = $i + 1;
        $j = $i - 1;

        $right = $positions[$k] - $positions[$i];
        $left = $positions[$i] - $positions[$j];

        while($right < $read_length)
        {
          $steps_right += 1;
          if($k == scalar(@positions) - 1){last;}
          $k += 1;
          $right = $positions[$k] - $positions[$i];  
        }

        while($left < $read_length)
        {
          $steps_left += 1;
          if($j == 0){last;}
          $j -= 1;
          $left = $positions[$i] - $positions[$j];  
        }

        $neighbors{$chr}{$positions[$i]} = $steps_right + $steps_left;
      }
    }
  }
  splice @positions;
}

open(OUT,"> $output") or die "Error writing to $output\n";

foreach $chr (keys %SNPhash)
{
  foreach $pos (keys %{$SNPhash{$chr}})
  {
    if($mappHash{$chr}{$pos})
    {
      #if($mappHash{$chr}{$pos} == 1 && $neighbors{$chr}{$pos} < $mis)
      if($mappHash{$chr}{$pos} == 1)
      {
        print OUT "$chr\t$pos\t$SNPhash{$chr}{$pos}[0]\t$SNPhash{$chr}{$pos}[1]\n";
      }
    }
  }
}
close OUT;
