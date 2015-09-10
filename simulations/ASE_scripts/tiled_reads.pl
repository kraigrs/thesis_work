#!/usr/bin/perl

########################################################################
# 
# 02/08/2012
#
# tiled_reads.pl
#
# Purpose: simulate SE or PE reads from a set from specified regions
#
# Notes: 10/27/2011: still working on how to represent the fact that
#                    genes will not be expressed in the same amount.
#                    
# 
# Input: type of read, length of sequence, error rate, both references, 
#        regions to generate reads from, levels of expression, and output location
#
# Fixes: 02/08/2012 --> 
#
# Output: FASTA files representing reads generated 
#
# Syntax: perl tiled_reads.pl single 50 0 ../../../Dmel/chromFa/chr2R.fa ../../../Dmel/chromFa/chr2R_alt.fa ../../sim_regs_chr2R.bed ../../../Simulations/tiled
#
########################################################################

use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use POSIX;

#srand(12345);

my$type = $ARGV[0];
my$length = $ARGV[1];
my$error = $ARGV[2]; # decimal 
my$ref = $ARGV[3];
my$refAlt = $ARGV[4];
my$regions = $ARGV[5]; # sim_regs_chr2L.bed 
my$out_dir = $ARGV[6];
 
my$refPos = ""; my$refNeg = ""; my$altPos = ""; my$altNeg = "";
our($chr,$pos,$code,$seq,$i,$j,$k,$reflength,$posNeg,$false,$base,$reg,$start,$stop,$begin,$end,$gene,$level,$readRef,$readCons);
our(@elements,@parts);
our(%regs,%expr);

print "Inputting $ref\n";

open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  #print "Here!\n";
  if(/^>chr(\S+)$/){$chr = $1;}
  else
  {
    $refPos = $refPos.uc($_); 
  }
}
close REF;
#exit;

$refNeg = &revComp($refPos);

#specify alternate genomes in both orientations
print "Inputting $refAlt\n";
open(REF,"$refAlt") or die "\nError opening $refAlt\n";
while(<REF>)    
{
  chomp;
  #print "Here!\n";
  if(/^>chr(\S+)$/){$chr = $1;}
  else
  {
    $altPos = $altPos.uc($_); 
  } 
}
close REF;

$altNeg = &revComp($altPos);

if(length($refPos) == length($altPos)){$reflength = length($refPos);}
else{die "The reference lengths don't match up!";}

# regions from which to generate reads
print "Inputting $regions\n";
open(REG,"$regions") or die "\nError opening $regions\n";
while(<REG>)    
{
  chomp;
  @elements = split("\t",$_); 
  $start = $elements[1];
  $stop = $elements[2];
  $gene = $elements[3];
  $reg = $start.",".$stop;
  #@parts = split("_",$elements[3]);
  #$gene = $parts[0];
  #$reg = $parts[1];
  #print "\nGene:$gene\tReg:$reg\tStart:$start\tStop:$stop";
  $regs{$gene}{$reg} = [$start,$stop];
}
close REG;

print "Generating reads for chr$chr\n\n";

# single end simulation
if($type eq "single")
{
  open(OUT,">$out_dir/chr$chr\_$type\_bp$length\_error$error\_tiled.fa") or die "\nError writing to $out_dir/chr$chr\_$type\_bp$length\_error$error\_tiled.fa\n";

  foreach $gene (keys %regs)
  { 
    #print "\nIn gene loop...";
    foreach $reg (keys %{$regs{$gene}})
    {
      if(($regs{$gene}{$reg}[1]-$regs{$gene}{$reg}[0]) >= $length)
      {
        # positive orientation

        @elements = ($regs{$gene}{$reg}[0]..($regs{$gene}{$reg}[1]-$length));
        foreach $i (@elements)
        {
          if($error == 0)
          {
            $readRef = substr $refPos, $i, $length;
            $readCons = substr $altPos, $i, $length;
          }
          else
          {
            # need to fix the error part
            $readRef = substr $refPos, $i, $length;
            $readCons = substr $altPos, $i, $length;

            $readRef = &baseError($readRef); 
            $readCons = &baseError($readCons);
          }

          $begin = $i; $end = $i+$length; #BED format

          print OUT "\>HWI_$chr\_$gene\_$begin\_$end\_refPos\n";
          print OUT "$readRef\n";
          print OUT "\>HWI_$chr\_$gene\_$begin\_$end\_altPos\n";
          print OUT "$readCons\n";
	}

        # negative orientation

        @elements = ( ($reflength-$regs{$gene}{$reg}[1])..($reflength-$regs{$gene}{$reg}[0]-$length));
        foreach $i (@elements) 
        {
          if($error == 0)
          {
            $readRef = substr $refNeg, $i, $length;
            $readCons = substr $altNeg, $i, $length;
          }
          else
          {
            # need to fix the error part
            $readRef = substr $refNeg, $i, $length;
            $readCons = substr $altNeg, $i, $length;

            $readRef = &baseError($readRef); 
            $readCons = &baseError($readCons);
          }
         
          $begin = $reflength-$i-$length; $end = $reflength-$i; #BED format. + orienation

          print OUT "\>HWI_$chr\_$gene\_$begin\_$end\_refNeg\n";
          print OUT "$readRef\n";
          print OUT "\>HWI_$chr\_$gene\_$begin\_$end\_altNeg\n";
          print OUT "$readCons\n";
	}
      }
    }
  }
  close OUT;
}

sub revComp
{
  my$val = shift;
  my$revComp = reverse $val;
  $revComp =~ tr/ACGT/TGCA/;
  return $revComp;
}

sub baseError
{
  my$base = shift;
  my$false = rand();
  my$type;
  if($false < $error)
  {
    $type = rand();
    if($type < 0.25){$base = "A";}
    elsif($type >= 0.25 && $type < 0.5){$base = "C";}
    elsif($type >= 0.5 && $type < 0.75){$base = "G";}
    else{$base = "T";}
  }
  return $base;
}

