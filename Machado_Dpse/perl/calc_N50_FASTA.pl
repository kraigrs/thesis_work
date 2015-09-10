#!/usr/bin/perl

########################################################################
# 
# 07/26/2013
#
# calc_N50_FASTA.pl
#
# Purpose: take in a reference genome and return definitions of N50
# 
# Input: FASTA 
#
# Output: N50's
#         
# Syntax: perl calc_N50_FASTA.pl <FASTA> 
#
########################################################################

use strict;
use warnings;

my$ref = $ARGV[0];

my$chr; my$seq; my$first; my$length; my$nuc; my$total; my$sum; my$frac; my$i; my$N50_ascend; my$N50_descend; my$N50_avg;
my@lengths; my@nuc; my@ascend; my@descend;
my%genome; my%ACGT;

# read in reference sequence

$seq="";
$first=1;

open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)\s*/)
  #if(/^\>(Unknown\S+)\s*/)
  {
    unless($first == 1){$genome{$chr} = $seq;} # the first time, don't save sequence in hash
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
$genome{$chr} = $seq; # this needs to be done since last entry will not be stored in hash

# calculate total genome length

foreach $chr (keys %genome)
{
  $length = length($genome{$chr});
  push @lengths, $length;
  $total += $length;
  print "$chr\t$length\n";

  #@nuc = split("",$genome{$chr});
  #foreach(@nuc){$ACGT{$_} += 1;}

  #for($i=0;$i<$length;$i++)
  #{
  #  $nuc = substr $genome{$chr},$i,1;
  #  $ACGT{$nuc} += 1;
  #}
}
print "Total: $total\n";
#foreach(keys %ACGT)
#{
#  print "Base: $_\tCount: $ACGT{$_}\n";
#}

# calculate N50 based on ascending order

@ascend = sort {$a <=> $b} @lengths;

$i = 0;
$sum = $ascend[$i];
$N50_ascend = $ascend[$i];

while($sum/$total <= 0.5)
{
  $N50_ascend = $ascend[$i];
  $i++;
  $sum += $ascend[$i];
} 

print "N50 ascending: $N50_ascend\n";

# calculate N50 based on descending order

@descend = sort {$b <=> $a} @lengths;

$i = 0;
$sum = $descend[$i];
$N50_descend = $descend[$i];

while($sum/$total <= 0.5)
{
  $N50_descend = $descend[$i];
  $i++;
  $sum += $descend[$i];
} 

print "N50 descending: $N50_descend\n";

if($N50_ascend != $N50_descend){$N50_avg = ($N50_ascend+$N50_descend)/2;}
else{$N50_avg = $N50_ascend;}

print "N50 average: $N50_avg\n";
