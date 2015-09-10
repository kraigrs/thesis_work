#!/usr/bin/perl

########################################################################
# 
# 01/09/2013 (will fill in later)
#
########################################################################

use strict;
use warnings;

my$bed1 = $ARGV[0];
my$bed2 = $ARGV[1];
my$genos = $ARGV[2];
my$out = $ARGV[3];
my$spec1 = $ARGV[4];
my$spec2 = $ARGV[5];

my$chr; my$start; my$stop; my$name; my$base; my$geno; my$pos; my$n;
my@elements; my@meta;
my%bed1; my%bed2; my%names; my%s1; my%s2; my%amb; my%err;

open(BED1,"$bed1") or die "\nError opening $bed1\n";
while(<BED1>)    
{
  chomp;
  @elements = split("\t",$_);
  if($elements[0] =~ /^chr(\S+)$/){$chr = $1;}
  $start = $elements[1];
  $stop = $elements[2];
  @meta = split(/\_/,$elements[3]);
  if(scalar(@meta) == 5){$name = $meta[0]."_".$meta[1]."_".$meta[2]."_".$meta[3];}
  else{$name = $meta[0]."_".$meta[1];}
  $base = $elements[4];
 
  $bed1{$chr}{$stop} = [$base,$name];
  $s1{$chr}{$stop}{$name} = 0;
  $amb{$chr}{$stop}{$name} = 0;
  $err{$chr}{$stop}{$name} = 0;
}
close BED1;

open(BED2,"$bed2") or die "\nError opening $bed2\n";
while(<BED2>)    
{
  chomp;
  @elements = split("\t",$_);
  if($elements[0] =~ /^chr(\S+)$/){$chr = $1;}
  $start = $elements[1];
  $stop = $elements[2];
  @meta = split(/\_/,$elements[3]);
  if(scalar(@meta) == 5){$name = $meta[0]."_".$meta[1]."_".$meta[2]."_".$meta[3];}
  else{$name = $meta[0]."_".$meta[1];}
  $base = $elements[4];
 
  $bed2{$chr}{$stop} = [$base,$name];
  $s2{$chr}{$stop}{$name} = 0;
  $amb{$chr}{$stop}{$name} = 0;
  $err{$chr}{$stop}{$name} = 0;
}
close BED2;

#foreach $chr (keys %bed1)
#{
#  foreach $pos (keys %{$bed1{$chr}})
#  {
#    $base = $bed1{$chr}{$pos}[0];
#    $name = $bed1{$chr}{$pos}[1];
#    print "$chr\t$pos\t$name\t$base\n";
#  }
#}
#exit;

print "\nReading DGRP...\n";

#$n = 0;
open(DGRP,"$genos") or die "\nError opening $genos\n";
while(<DGRP>)
{
  chomp;
  @elements = split(/,/,$_);
  $chr = $elements[0];
  $pos = $elements[1];
  $geno = $elements[3];
  #print "$chr\t$pos\t$geno\n";  

  #$n += 1;
  #if(int($n/1000000) == $n/1000000){print "On line $n\n";}

  if($bed1{$chr}{$pos} && $bed2{$chr}{$pos} && $bed1{$chr}{$pos}[1] eq $bed2{$chr}{$pos}[1])
  {
    $name = $bed1{$chr}{$pos}[1];
    $names{$chr}{$pos} = $name; 

    if($geno =~ /[^ACGTacgt]/)
    {
      $amb{$chr}{$pos}{$name} += 1;
    }
    elsif($geno eq $bed1{$chr}{$pos}[0])
    {
      $s1{$chr}{$pos}{$name} += 1;
    }
    elsif($geno eq $bed2{$chr}{$pos}[0])
    {
      $s2{$chr}{$pos}{$name} += 1;
    }
    else
    {
      $err{$chr}{$pos}{$name} += 1;
    }
  }
}
close DGRP;

#print "\nDone!\n";

open(OUT,"\> $out") or die "Error opening $out!";
print OUT "chr\tposition\tlocus\t$spec1\t$spec2\tambiguous\terror\n";

foreach $chr (keys %names)
{
  foreach $pos (keys %{$names{$chr}})
  {
    $name = $names{$chr}{$pos};
    print OUT "chr$chr\t$pos\t$name\t$s1{$chr}{$pos}{$name}\t$s2{$chr}{$pos}{$name}\t$amb{$chr}{$pos}{$name}\t$err{$chr}{$pos}{$name}\n";
  }
}
close OUT;
