#!/usr/bin/perl

########################################################################
# 
# 01/09/2013 (will fill in later)
#
# Example: perl make_SNP_BED.pl Aldh/regions_zhr_Aldh.bed Aldh/regions_z30_Aldh.bed Aldh/regions_zhr_z30_Aldh.seqs_fsa.SNPs.txt Aldh/regions_zhr_Aldh.SNPs.bed Aldh/regions_z30_Aldh.SNPs.bed
#
########################################################################

use strict;
use warnings;

my$regions1 = $ARGV[0];
my$regions2 = $ARGV[1];
my$SNPs = $ARGV[2];
my$bed1 = $ARGV[3];
my$bed2 = $ARGV[4];

my$chr; my$start; my$stop; my$name; my$pos1; my$base1; my$pos2; my$base2; my$link;
my@elements;
my%regs1; my%regs2; my%polys; my%names;

open(REG1,"$regions1") or die "\nError opening $regions1\n";
while(<REG1>)    
{
  chomp;
  @elements = split("\t",$_);
  $chr = $elements[0];
  $start = $elements[1];
  $stop = $elements[2];
  $name = $elements[3];
 
  $regs1{$name} = [$chr,$start,$stop];
  $names{$name} = 1;
}
close REG1;

open(REG2,"$regions2") or die "\nError opening $regions2\n";
while(<REG2>)    
{
  chomp;
  @elements = split("\t",$_);
  $chr = $elements[0];
  $start = $elements[1];
  $stop = $elements[2];
  $name = $elements[3];
 
  $regs2{$name} = [$chr,$start,$stop];
  $names{$name} = 1;
}
close REG2;

#print "\n";
#foreach $name (keys %regs)
#{
#  $chr = $regs{$name}[0];
#  $start = $regs{$name}[1];
#  $stop = $regs{$name}[2];
#  print "$name\t$chr\t$start\t$stop\n";
#}
#print "\n";

open(SNP,"$SNPs") or die "\nError opening $SNPs\n";
while(<SNP>)    
{
  chomp;
  @elements = split("\t",$_);
  $name = $elements[0];
  $pos1 = $elements[1];
  $pos2 = $elements[2];
  $link = $elements[3];
  $base1 = $elements[4];
  $base2 = $elements[5];

  $polys{$name}{$link} = [$pos1,$base1,$pos2,$base2];
  $names{$name} = 1;
}
close SNP;

open(BED1,"\> $bed1") or die "Error opening $bed1!";
open(BED2,"\> $bed2") or die "Error opening $bed2!";

foreach $name (keys %names)
{
  if($regs1{$name} && $regs2{$name} && $polys{$name})
  {
    #print "Here";
    foreach $link (keys %{$polys{$name}})
    {
      $chr = $regs1{$name}[0];
      $pos1 = $polys{$name}{$link}[0];
      $base1 = $polys{$name}{$link}[1];
      $start = $regs1{$name}[1]+$pos1-1;
      $stop = $start+1;
      print BED1 "$chr\t$start\t$stop\t$name\_$link\t$base1\n";

      $chr = $regs2{$name}[0];
      $pos2 = $polys{$name}{$link}[2];
      $base2 = $polys{$name}{$link}[3];
      $start = $regs2{$name}[1]+$pos2-1;
      $stop = $start+1;
      print BED2 "$chr\t$start\t$stop\t$name\_$link\t$base2\n";

    }
  }
}
close BED1; close BED2;
