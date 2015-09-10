#!/usr/bin/perl

###################################################################################
#
# 07/12/2011
#
# hetSites.pl
#
# Purpose: get the commons heterozygous sites for zhr and z30
# 
# Input: 4 het files (BED) in dm3 space
#
# Output: common set of heterozygous SNPs and INDELs 
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
  my$SNP1 = $ARGV[0];
  my$SNP2 = $ARGV[1];
  my$INDEL1 = $ARGV[2];
  my$INDEL2 = $ARGV[3];
  my$SNPout = $ARGV[4];
  my$INDELout = $ARGV[5];

  our($chr,$start,$stop);
  our(@elements);
  our(%hetSNPs,%hetINDELs);

  open(SNP1,"$SNP1") or die "Can't open $SNP1 for reading!";
  while(<SNP1>) 
  {
    chomp;
    @elements = split(/\s+/,$_);
    $chr = $elements[0];
    $start = $elements[1];
    $stop = $elements[2]; 
    $hetSNPs{$chr}{$start} = $stop;
  }
  close SNP1;

  open(SNP2,"$SNP2") or die "Can't open $SNP2 for reading!";
  while(<SNP2>) 
  {
    chomp;
    @elements = split(/\s+/,$_);
    $chr = $elements[0];
    $start = $elements[1];
    $stop = $elements[2];
    $hetSNPs{$chr}{$start} = $stop;
  }
  close SNP2;

  open(INDEL1,"$INDEL1") or die "Can't open $INDEL1 for reading!";
  while(<INDEL1>) 
  {
    chomp;
    @elements = split(/\s+/,$_);
    $chr = $elements[0];
    $start = $elements[1];
    $stop = $elements[2];
    $hetINDELs{$chr}{$start} = $stop;
  }
  close INDEL1;
  
  open(INDEL2,"$INDEL2") or die "Can't open $INDEL2 for reading!";
  while(<INDEL2>) 
  {
    chomp;
    @elements = split(/\s+/,$_);
    $chr = $elements[0];
    $start = $elements[1];
    $stop = $elements[2];
    $hetINDELs{$chr}{$start} = $stop;
  }
  close INDEL2;

  open(OUT,">$SNPout") or die "Error writing to $SNPout\n";
  foreach $chr (keys %hetSNPs)
  {
    foreach $start (keys %{$hetSNPs{$chr}})
    {
      print OUT "$chr\t$start\t$hetSNPs{$chr}{$start}\n";
    }  
  }
  close OUT;

  open(OUT,">$INDELout") or die "Error writing to $INDELout\n";
  foreach $chr (keys %hetINDELs)
  {
    foreach $start (keys %{$hetINDELs{$chr}})
    {
      print OUT "$chr\t$start\t$hetINDELs{$chr}{$start}\n";
    }  
  }
  close OUT;

}
