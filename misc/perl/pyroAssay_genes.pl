#!/usr/bin/perl

###################################################################################
#
# 11/18/2010
#
# pyroAssay_genes.pl
#
# Purpose: create sequences to be used in a pyrosequencing assay (i.e. from zhr and z30) 
# 
# Input: list of loci of interest in dm3, coordinates mapped to zhr and z30
#
# Output: for each dm3 locus, the correponding zhr and z30 sequences 
#
# Usage: perl pyroAssay.pl <list> <constExonsS1> <constExonsS2> <ref1> <ref2> <output>  
#
#        <list>         ==> a list of genes
#        <constExonsS1> ==> a file containing a list of constitutive exons from species 1
#        <constExonsS2> ==> a file containing a list of constitutive exons from species 2
#        <ref1>         ==> fasta reference from species 1
#        <ref2>         ==> fasta reference from species 2
#        <output>       ==> the output file to write converted reads to
#
# Example: perl pyroAssay_genes.pl ../mel_mel_data/imprinting_list2_genes_062711.txt ../McManus/constExons.zhr.bed ../McManus/constExons.z30.bed ../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_AMBIG_all.fa ../mel_mel_data/Resequencing/resequencing-assembly/z30/z30_AMBIG_all.fa ../mel_mel_data/imprinting_list2_pyroSeqs_062711.txt
#
###################################################################################

use strict;
use warnings;

main ();
sub main 
{
  my$list = $ARGV[0];  
  my$constExonsS1 = $ARGV[1]; 
  my$constExonsS2 = $ARGV[2];
  my$ref1 = $ARGV[3];
  my$ref2 = $ARGV[4];
  my$output = $ARGV[5];

  my$chr; my$segment; my$seqFound = 0;
  my%ref1; my%ref2;

  open(REF1,"$ref1") or die "Can't open $ref1 for reading!\n";
  while(<REF1>) 
  {
    chomp;
    if(/^>(\S+)/)
    {
      if($seqFound == 1)
      {
        $ref1{$chr} = $segment;
        $seqFound = 0;
        undef $segment;
      }
      $chr = $1;
    }
    else
    {
      $seqFound = 1;
      $segment = $segment . $_;
      $ref1{$chr} = $segment;
    }
  }
  close REF1;
  undef $segment; $seqFound = 0;

  open(REF2,"$ref2") or die "Can't open $ref2 for reading!\n";
  while(<REF2>) 
  {
    chomp;
    if(/^>(\S+)/)
    {
      if($seqFound == 1)
      {
        $ref2{$chr} = $segment;
        $seqFound = 0;
        undef $segment;
      }
      $chr = $1;
    }
    else
    {
      $seqFound = 1;
      $segment = $segment . $_;
      $ref2{$chr} = $segment;
    }
  }
  close REF2;
  undef $segment; $seqFound = 0;

  my@elements; my@locus;
  my%exons1; my%exons2;

  open(EXONS1,"$constExonsS1") or die "Can't open $constExonsS1 for reading!\n";
  while(<EXONS1>) 
  {
    chomp;
    unless (/^#track.+/)
    {
      @elements = split(/\s+/,$_); 
      @locus = split(/_/,$elements[3]);
      $exons1{$locus[0]}{$locus[1]} = [$elements[0],$elements[1],$elements[2]];
      #exons1{gene}{exon} = [location, start, stop]
    }
  }
  close EXONS1;

  open(EXONS2,"$constExonsS2") or die "Can't open $constExonsS2 for reading!\n";
  while(<EXONS2>) 
  {
    chomp;
    unless (/^#track.+/)
    {
      @elements = split(/\s+/,$_); 
      @locus = split(/_/,$elements[3]);
      $exons2{$locus[0]}{$locus[1]} = [$elements[0],$elements[1],$elements[2]];
      #exons2{gene}{exon} = [location, start, stop]
    }
  }
  close EXONS2;
  #print "\nDone making hashes!\n";

  my$start; my$stop; my$i; my$k1; my$k2;
  my@seq;

  open(OUT,">$output") or die "Error writing to $output!\n";

  open (LIST,"$list") or die "Can't open $list for reading!\n";
  while (<LIST>) 
  {
    chomp;
    $k1 = $_;
    for $k2 (keys %{$exons1{$k1}})
    {
      if($exons2{$k1}{$k2})
      {
        print OUT "\>$k1\_$k2\,$exons1{$k1}{$k2}[0]\_$exons1{$k1}{$k2}[1]\,$exons1{$k1}{$k2}[2]\n";

        $start = $exons1{$k1}{$k2}[1];
        $stop = $exons1{$k1}{$k2}[2];
        @seq = split("",$ref1{$exons1{$k1}{$k2}[0]});

        for ($i = $start; $i < $stop; $i++) 
        {
          print OUT "$seq[$i]";
        }

        print OUT "\n\>$k1\_$k2\,$exons2{$k1}{$k2}[0]\_$exons2{$k1}{$k2}[1]\,$exons2{$k1}{$k2}[2]\n";

        $start = $exons2{$k1}{$k2}[1];
        $stop = $exons2{$k1}{$k2}[2];
        @seq = split("",$ref2{$exons2{$k1}{$k2}[0]});

        for ($i = $start; $i < $stop; $i++) 
        {
 	  print OUT "$seq[$i]";
        }
        print OUT "\n";
      }
    }
  }
  close LIST;
  close OUT;
}
