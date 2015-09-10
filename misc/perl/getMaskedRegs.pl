#!/usr/bin/perl

###################################################################################
#
# 08/13/2010
#
# getMaskedRegs.pl
#
# Purpose: using individual parents aligned to both genomes, create a gap file in dm3
#          representing regions that misclassify reads  
# 
# Input: the classified reads from classify.pl, probably called *readout*, and the lifted (and converted)
#        .bed from the alignment files  
#
# Output: a list of regions that are "gaps", that is, regions that should be masked due to incorrect alignment
#
# E.g. perl separate.pl <species1> <species2> <readout> 
#
#      <species1>  ==> indicate name of correct species (i.e. Dmel, Dsec, etc.)
#      <species2>  ==> indicate name of other species aligned to
#      <readoutS1> ==> readout file from classify.pl
#      <readoutS2> ==> readout file from classify.pl
#      <s1mate1al2>   ==> a lifted .bed file from the 1st mate of the 1st species aligned to 2nd genome
#      <s1mate2al2>   ==> a lifted .bed file from the 2nd mate of the 1st species aligned to 2nd genome
#      <s2mate1al1>   ==> a lifted .bed file from the 1st mate of the 2nd species aligned to 1st genome
#      <s2mate2al1>   ==> a lifted .bed file from the 2nd mate of the 2nd species aligned to 1st genome
#      <output>    ==> prefix for output files
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
  my$species1 = $ARGV[0]; # zhr
  my$species2 = $ARGV[1]; # z30
  my$readoutS1 = $ARGV[2];
  my$readoutS2 = $ARGV[3];
  my$s1mate1al2 = $ARGV[4];
  my$s1mate2al2 = $ARGV[5];
  my$s2mate1al1 = $ARGV[6];
  my$s2mate2al1 = $ARGV[7];
  my$output = $ARGV[8]; 

  my$read; my$mate1S1; my$mate1S2; my$mate2S1; my$mate2S2; my$class;
  my@elements;
  my%gapS1mate1al2; my%gapS1mate2al2; my%gapS2mate1al1; my%gapS2mate2al1;
   
  # open readoutS1 file
  open(READS1,"$readoutS1") or die "Can't open $readoutS1 for reading!";
  LINE:while(<READS1>) 
  {
    chomp;  
    if($_ =~ /^read/){next LINE;}
    elsif($_ =~ /^HWI.+/)
    { 
      #print "$_\n";
      @elements = split(/\s/,$_);
      $read = $elements[0];
      $mate1S1 = $elements[1];
      $mate1S2 = $elements[2];
      $mate2S1 = $elements[3];
      $mate2S2 = $elements[4];
      $class = $elements[5];
      #print "$class\n";
      if($class eq $species2)
      {
        if($mate1S1 ne $mate1S2 && $mate2S1 ne $mate2S2){$gapS1mate1al2{$read} = 1;$gapS1mate2al2{$read} = 1;}
        elsif($mate1S1 ne $mate1S2){$gapS1mate1al2{$read} = 1;}
        elsif($mate2S1 ne $mate2S2){$gapS1mate2al2{$read} = 1;}
      }  
    }
  }
  close READS1;

  # open readoutS2 file
  open(READS2,"$readoutS2") or die "Can't open $readoutS2 for reading!";
  LINE:while(<READS2>) 
  {
    chomp;  
    if($_ =~ /^read/){next LINE;}
    elsif($_ =~ /^HWI.+/)
    { 
      #print "$_\n";
      @elements = split(/\s/,$_);
      $read = $elements[0];
      $mate1S1 = $elements[1];
      $mate1S2 = $elements[2];
      $mate2S1 = $elements[3];
      $mate2S2 = $elements[4];
      $class = $elements[5];
      #print "$class\n";
      if($class eq $species1)
      {
        if($mate1S1 ne $mate1S2 && $mate2S1 ne $mate2S2){$gapS2mate1al1{$read} = 1;$gapS2mate2al1{$read} = 1;}
        elsif($mate1S1 ne $mate1S2){$gapS2mate1al1{$read} = 1;}
        elsif($mate2S1 ne $mate2S2){$gapS2mate2al1{$read} = 1;}
      }  
    }
  }
  close READS2;

  print "\nDone making reads hashes!\n";

  my$chr; my$start; my$stop;

  open(GAPS,">$output") or die "Error writing to $output\n";

  open(S1MATE1AL2,"$s1mate1al2") or die "Can't open $s1mate1al2 for reading!";
  while(<S1MATE1AL2>) 
  {
    chomp;
    @elements = split(/\s/,$_);
    $chr = $elements[0];
    $start = $elements[1];
    $stop = $elements[2];
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}

    if($gapS1mate1al2{$read}){print GAPS "$chr\t$start\t$stop\n";}
  }
  close S1MATE1AL2;

  open(S1MATE2AL2,"$s1mate2al2") or die "Can't open $s1mate2al2 for reading!";
  while(<S1MATE2AL2>) 
  {
    chomp;
    @elements = split(/\s/,$_);
    $chr = $elements[0];
    $start = $elements[1];
    $stop = $elements[2];
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}

    if($gapS1mate2al2{$read}){print GAPS "$chr\t$start\t$stop\n";}
  }
  close S1MATE2AL2;
    
  open(S2MATE1AL1,"$s2mate1al1") or die "Can't open $s2mate1al1 for reading!";
  while(<S2MATE1AL1>) 
  {
    chomp;
    @elements = split(/\s/,$_);
    $chr = $elements[0];
    $start = $elements[1];
    $stop = $elements[2];
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}

    if($gapS2mate1al1{$read}){print GAPS "$chr\t$start\t$stop\n";}
  }
  close S2MATE1AL1;

  open(S2MATE2AL1,"$s2mate1al1") or die "Can't open $s2mate1al1 for reading!";
  while(<S2MATE2AL1>) 
  {
    chomp;
    @elements = split(/\s/,$_);
    $chr = $elements[0];
    $start = $elements[1];
    $stop = $elements[2];
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}

    if($gapS2mate2al1{$read}){print GAPS "$chr\t$start\t$stop\n";}
  }
  close S2MATE2AL1; 
  close GAPS;
} 
 
