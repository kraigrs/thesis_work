#!/usr/bin/perl

###################################################################################
#
# 11/18/2010
#
# pyroAssay_exons.pl
#
# Purpose: create sequences to be used in a pyrosequencing assay (i.e. from zhr and z30) 
# 
# Input: list of loci of interest in dm3, coordinates mapped to zhr and z30
#
# Output: for each dm3 locus, the correponding zhr and z30 sequences 
#
# Usage: perl pyroAssay.pl <list> <constExonsS1> <constExonsS2> <ref1> <ref2> <output>  
#
#        <list>         ==> a list of constitutive exons in the following format: name_start,stop (in parent species space)
#        <constExonsS1> ==> a file containing a list of constitutive exons from species 1
#        <constExonsS2> ==> a file containing a list of constitutive exons from species 2
#        <ref1>         ==> fasta reference from species 1
#        <ref2>         ==> fasta reference from species 2
#        <output>       ==> the output file to write converted reads to
#
# Example: perl pyroAssay.pl ../mel_mel_data/imprinting_exons_for_pyro_062711.txt ../McManus/constExons.zhr.bed ../McManus/constExons.z30.bed ../mel_mel_data/Resequencing/resequencing-assembly/zhr/zhr_AMBIG_all.fa ../mel_mel_data/Resequencing/resequencing-assembly/z30/z30_AMBIG_all.fa ../mel_mel_data/imprinting_pyro_seqs_062711.txt
#
###################################################################################

use strict;
use warnings;

main ();
sub main 
{
  my$list = $ARGV[0];  
  my$ref1 = $ARGV[1];
  my$constExonsS1 = $ARGV[2]; 
  my$ref2 = $ARGV[3];
  my$constExonsS2 = $ARGV[4];
  my$output = $ARGV[5];

  my$chr; my$first = 1; my$seq = ""; my$length;
  my%ref1; my%ref2;

  open(REF1,"$ref1") or die "\nError opening $ref1\n";
  while(<REF1>)    
  {
    chomp;
    if(/^\>(\S+)$/)
    {
      unless($first == 1){$ref1{$chr} = $seq;} # the first time, don't save sequence in hash
      $chr = $1;
      $seq = "";
      $first = 0;
    }
    else
    {
      $seq = $seq.uc($_);
    }
  }
  close REF1;
  $ref1{$chr} = $seq; # this needs to be done since last entry will not be stored in hash

  $first = 1; $seq = "";
  open(REF2,"$ref2") or die "\nError opening $ref2\n";
  while(<REF2>)    
  {
    chomp;
    if(/^\>(\S+)$/)
    {
      unless($first == 1){$ref2{$chr} = $seq;} # the first time, don't save sequence in hash
      $chr = $1;
      $seq = "";
      $first = 0;
    }
    else
    {
      $seq = $seq.uc($_);
    }
  }
  close REF2;
  $ref2{$chr} = $seq; # this needs to be done since last entry will not be stored in hash

  #foreach $chr (keys %ref1)
  #{
  #  $length = length($ref1{$chr});
  #  print "Chromosome: $chr\tLength: $length\n";
  #}

  #foreach $chr (keys %ref2)
  #{
  #  $length = length($ref2{$chr});
  #  print "Chromosome: $chr\tLength: $length\n";
  #}
  #exit;

  my@elements;
  my%exons1; my%exons2;

  open(EXONS1,"$constExonsS1") or die "Can't open $constExonsS1 for reading!\n";
  while(<EXONS1>) 
  {
    chomp;
    unless (/^#track.+/)
    {
      @elements = split(/\s+/,$_); 
      $exons1{$elements[3]} = [$elements[0],$elements[1],$elements[2]];
      #exons1{chrm} = [location, start, stop]
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
      $exons2{$elements[3]} = [$elements[0],$elements[1],$elements[2]];
      #exons2{chrm} = [location, start, stop]
    }
  }
  close EXONS2;
  #print "\nDone making hashes!\n";

  my$start; my$stop; #my$i;
  #my@seq;

  open(OUT,">$output") or die "Error writing to $output!\n";

  open(LIST,"$list") or die "Can't open $list for reading!\n";
  while (<LIST>) 
  {
    chomp;
    if($exons1{$_} && $exons2{$_})
    {
      #print "Here!";
      print OUT "\>$_\:$exons1{$_}[0]\_$exons1{$_}[1]\_$exons1{$_}[2]\n";

      $chr = $exons1{$_}[0];
      $start = $exons1{$_}[1];
      $stop = $exons1{$_}[2];
      $length = $stop - $start;      

      $seq = substr $ref1{$chr},$start,$length;
      print OUT "$seq";

      #@seq = split("",$ref1{$exons1{$_}[0]});
      #for ($i = $start; $i < $stop; $i++) 
      #{
      #  print OUT "$seq[$i]";
      #}

      print OUT "\n\>$_\:$exons2{$_}[0]\_$exons2{$_}[1]\_$exons2{$_}[2]\n";

      $chr = $exons2{$_}[0];
      $start = $exons2{$_}[1];
      $stop = $exons2{$_}[2];
      $length = $stop - $start;      

      $seq = substr $ref2{$chr},$start,$length;
      print OUT "$seq";

      #@seq = split("",$ref2{$exons2{$_}[0]});
      #for ($i = $start; $i < $stop; $i++) 
      #{
      #	 print OUT "$seq[$i]";
      #}
      print OUT "\n";
    }
  }
  close LIST;
  close OUT;
}
