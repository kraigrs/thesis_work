#!/usr/bin/perl

########################################################################
# 
# 08/01/2013
#
# get_orthologs.pl
#
# Purpose: curate list of orthologs between Dmel and Dpse
# 
# Input: Dmel and Dpse gene annotations 
#
# Output: one gene per line with Dmel and Dpse information
#         
# Syntax: perl get_orthologs.pl <dmel> <dpse>
#
########################################################################

use strict;
use warnings;

my$dmel = $ARGV[0];
my$dpse = $ARGV[1];

my$line; my$gene; my$FBid; my$chr; my$start; my$stop; my$strand;
my@elements;
my%mel2pse; my%mel_info; my%pse; my%pse_info;

open(DMEL,"$dmel") or die "\nError opening $dmel\n";
while(<DMEL>)    
{
  chomp;
  $line = $_;

  if($line =~ /^FBgn/)
  {
    @elements = split("\t",$line);
    if($elements[3] =~ /Dpse\\(GA\d+)/)
    {
      $gene = $1;
      $FBid = $elements[0];
      #$chr = "chr".$elements[5];
      $chr = $elements[5];
      $start = $elements[6]-1; #BED-style coordinates
      $stop = $elements[7];
      if($elements[8] == 1){$strand = "+";}
      else{$strand = "-";}

      $mel2pse{$gene} = $FBid;
      $mel_info{$FBid} = [$chr,$start,$stop,$strand];

      #print "$FBid\t$gene\t$chr\t$start\t$stop\t$strand\n";
    }
  }
}
close DMEL;

open(DPSE,"$dpse") or die "\nError opening $dpse\n";
while(<DPSE>)    
{
  chomp;
  $line = $_;

  if($line =~ /^FBgn/)
  {
    @elements = split("\t",$line);

    $FBid = $elements[0];
    $gene = $elements[2];
    #if($elements[5] =~ /Unknown/){$chr = $elements[5];}
    #else{$chr = "chr".$elements[5];}
    $chr = $elements[5];
    $start = $elements[6]-1; #BED-style coordinates
    $stop = $elements[7];
    if($elements[8] == 1){$strand = "+";}
    else{$strand = "-";}

    $pse{$gene} = $FBid;
    $pse_info{$FBid} = [$chr,$start,$stop,$strand];

    #print "$FBid\t$gene\t$chr\t$start\t$stop\t$strand\n";
  }
}
close DPSE;

print "dpse_FBid\tchr\tstart\tstop\tstrand\tdmel_FBid\tchr\tstart\tstop\tstrand\n";

foreach $gene (keys %pse)
{
  if($mel2pse{$gene})
  {
    print "$pse{$gene}\t";
    print "$pse_info{$pse{$gene}}[0]\t";
    print "$pse_info{$pse{$gene}}[1]\t";
    print "$pse_info{$pse{$gene}}[2]\t";
    print "$pse_info{$pse{$gene}}[3]\t";
    print "$mel2pse{$gene}\t";
    print "$mel_info{$mel2pse{$gene}}[0]\t";
    print "$mel_info{$mel2pse{$gene}}[1]\t";
    print "$mel_info{$mel2pse{$gene}}[2]\t";
    print "$mel_info{$mel2pse{$gene}}[3]\n";
  }
}
