#!/usr/bin/perl

########################################################################
# 
# 02/05/2014
#
# regulator_properties.pl
#
# Purpose: summarize the inferred gene regulatory network from Marbach et al.
# 
# Input: list of regulators, target genes, and score connecting them 
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
my$regdiv = $ARGV[2];

my$regulator; my$target; my$score; my$n; my$sum; my$avg; my$tot; my$counter;
my$FBgn; my$CG; my$class; my$trans; my$gene; my$sum_CG; my$avg_CG; my$prop_trans;
my@elements;
my%network; my%support; my%FBgn2CG; my%divergence;

open(CONN,"$connections") or die "\nError opening $connections\n";
while(<CONN>)
{
  chomp;
  @elements = split(/\s+/,$_);
  $regulator = $elements[0];
  $target = $elements[1];
  $score = $elements[2];
  #print "$regulator\t$target\t$score\n";

  push( @{ $network{$regulator} }, $target);
  push( @{ $support{$regulator} }, $score);
}
close CONN;

open(ANNO,"$annotations") or die "\nError opening $annotations\n";
while(<ANNO>)
{
  chomp;
  unless(/^\#/)
  {
    @elements = split(/\s+/,$_);
    $FBgn = $elements[0];
    $CG = $elements[2];

    $FBgn2CG{$FBgn} = $CG;
  }
}
close ANNO;

open(REG,"$regdiv") or die "\nError opening $regdiv\n";
while(<REG>)
{
  chomp;
  unless(/^gene/)
  {
    @elements = split(/\s+/,$_);
    $gene = $elements[0];

    if($elements[9] == 1){$class = "cis";}
    elsif($elements[10] == 1){$class = "trans";}
    elsif($elements[11] == 1){$class = "conserved";}
    elsif($elements[12] == 1){$class = "compensatory";}
    elsif($elements[13] == 1){$class = "ambiguous";}
    elsif($elements[14] == 1){$class = "cis+trans";}
    elsif($elements[15] == 1){$class = "cisXtrans";}
    else{$class = "none";}

    $divergence{$gene} = $class;
  }
}
close REG;

print "regulator\ttargets\tavg_score\tFB_targets\treg_div\ttrans\n";

foreach $regulator (keys %network)
{

  if($FBgn2CG{$regulator} && $divergence{$FBgn2CG{$regulator}})
  #if($FBgn2CG{$regulator})
  {
    # all targets
    $sum = 0;
    $tot = scalar( @{ $network{$regulator} } );
    foreach( @{ $support{$regulator} } ){$sum += $_;}
    $avg = $sum/$tot;

    # targets in expression data
    $n = 0;
    $trans = 0;
    foreach( @{ $network{$regulator} } )
    {
      if($FBgn2CG{$_} && $divergence{$FBgn2CG{$_}})
      {
        $n += 1;
        $gene = $FBgn2CG{$_};

        if($divergence{$gene} eq "trans"){$trans += 1;}
        elsif($divergence{$gene} eq "cis+trans"){$trans += 1;}
        elsif($divergence{$gene} eq "cisXtrans"){$trans += 1;}
        elsif($divergence{$gene} eq "compensatory"){$trans += 1;}
      }
    }
    
    if($n > 0){$prop_trans = $trans/$n;}
    else{$prop_trans = 0;}

    if($divergence{$FBgn2CG{$regulator}}){print "$FBgn2CG{$regulator}\t$tot\t$avg\t$n\t$divergence{$FBgn2CG{$regulator}}\t$trans\n";}
    else{print "$FBgn2CG{$regulator}\t$tot\t$avg\t$n\tnone\t$trans\n";}
  }
}
