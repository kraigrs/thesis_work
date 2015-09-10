#!/usr/bin/perl

########################################################################
# 
# 01/13/2013 (will fill in later)
#
########################################################################

use strict;
use warnings;

my$spec1 = $ARGV[0];
my$ref   = $ARGV[1]; #z30

my$spec2 = $ARGV[2];
my$file1 = $ARGV[3]; #zhr

my$spec3 = $ARGV[4];
my$file2 = $ARGV[5]; #sim

my$spec4 = $ARGV[6];
my$file3 = $ARGV[7]; #sec

my$out   = $ARGV[8];

my$name; my$base; my$locus; my$pos;
my$major; my$a1; my$a2; my$a3;
my@elements;
my%ref; my%s1; my%s2; my%s3; my%loci;

#z30
open(BED,"$ref") or die "\nError opening $ref\n";
while(<BED>)    
{
  chomp;
  @elements = split("\t",$_);
  
  $name = $elements[3];
  $base = $elements[4];
 
  $ref{$name} = $base;
  $loci{$name} = 1;
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
 
  $s1{$name} = $base;
  $loci{$name} = 1;
}
close BED;

#sim
open(BED,"$file2") or die "\nError opening $file2\n";
while(<BED>)    
{
  chomp;
  @elements = split("\t",$_);
  
  $name = $elements[3];
  $base = $elements[4];
 
  $s2{$name} = $base;
  $loci{$name} = 1;
}
close BED;

#sec
open(BED,"$file3") or die "\nError opening $file3\n";
while(<BED>)    
{
  chomp;
  @elements = split("\t",$_);
  
  $name = $elements[3];
  $base = $elements[4];
 
  $s3{$name} = $base;
  $loci{$name} = 1;
}
close BED;

# compare alleles

open(OUT,"> $out") or die "Error writing to $out\n";
print OUT "region\tpos\t$spec1\t$spec2\t$spec3\t$spec4\n";

foreach $name (keys %loci)
{
  if($ref{$name} && $s1{$name} && $s2{$name} && $s3{$name})
  {
    if($name =~ /^(\S+)\_(\d+)$/){$locus = $1; $pos = $2;}

    $major = $ref{$name};
    $a1 = $s1{$name};
    $a2 = $s2{$name};
    $a3 = $s3{$name};

    if($major ne $a1 && $a1 eq $a2 && $a2 eq $a3)
    {
      print OUT "$locus\t$pos\t$ref{$name}\t$s1{$name}\t$s2{$name}\t$s3{$name}\n";
    }
  }
}
