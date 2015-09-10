#!/usr/bin/perl

########################################################################
# 
# 04/17/2014
#
# chromosomal_autoregulation.pl
#
# Purpose: for regulators on different chromosomes, how many of their targets are on same chromosome
# 
# Input: connections (regulator --> target), gene annotations, subset of genes to consider 
#
# Output: for each regulator --> number of targets and average score for all targets
#         
# Syntax: perl regulator_properties.pl <list> 
#
########################################################################

use strict;
use warnings;

my$connections = $ARGV[0];
my$annotations = $ARGV[1];
#my$subset = $ARGV[2];

my$regulator; my$target; my$score; my$n; 
my$FBgn; my$CG; my$gene; my$chr; my$autoreg;
my@elements;
my%network; my%support; my%FBgn2CG; my%set; my%chroms;

open(CONN,"$connections") or die "\nError opening $connections\n";
while(<CONN>)
{
  chomp;
  @elements = split(/\s+/,$_);
  $regulator = $elements[0];
  $target = $elements[2];
  #$score = $elements[2];
  #print "$regulator\t$target\t$score\n";

  push( @{ $network{$regulator} }, $target);
  #push( @{ $support{$regulator} }, $score);
}
close CONN;

open(ANNO,"$annotations") or die "\nError opening $annotations\n";
while(<ANNO>)
{
  chomp;
  #unless(/^\#/)
  unless(/^\Ensembl/)
  {
    @elements = split(/\s+/,$_);
    $FBgn = $elements[2];
    #$CG = $elements[2];
    $chr = $elements[1];

    #$FBgn2CG{$FBgn} = $CG;
    $chroms{$FBgn} = $chr;
  }
}
close ANNO;

#open(SET,"$subset") or die "\nError opening $subset\n";
#while(<SET>)
#{
#  chomp;
#  unless(/^gene/)
#  {
#    $gene = $_;
#    #print "\{\{\{$gene\}\}\}\n";
#    $set{$gene} = 1;
#  }
#}
#close SET;

print "regulator\tchromosome\ttargets\tautoreg\n";

foreach $regulator (keys %network)
{
  #print "$regulator\t$FBgn2CG{$regulator}\n";
  #print "$set{$FBgn2CG{$regulator}}\n";

  #if($FBgn2CG{$regulator} && $set{$FBgn2CG{$regulator}})
  if($chroms{$regulator})
  {
    #print "here\n";
    $n = 0;
    $autoreg = 0;

    foreach $target ( @{ $network{$regulator} } )
    {
      #if($FBgn2CG{$target} && $set{$FBgn2CG{$target}})
      if($chroms{$target})
      {
        $n += 1;

        if($chroms{$regulator} eq $chroms{$target}){$autoreg += 1;}
      }
    }

    #print "$FBgn2CG{$regulator}\t$chroms{$regulator}\t$n\t$autoreg\n";
    print "$regulator\t$chroms{$regulator}\t$n\t$autoreg\n";
  }
}
