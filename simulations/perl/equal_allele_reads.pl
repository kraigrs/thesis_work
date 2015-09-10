#!/usr/bin/perl

########################################################################
# 
# 09/23/2011
#
# equal_allele_reads.pl
#
# Purpose: simulate SE or PE reads from a set from specified regions
#
# Notes: 10/27/2011: this script comes closer to reality by selecting
#                    a set of reads based on expression values from the DGRP
# 
# Input: type of read, length of sequence, error rate, both references, 
#        regions to generate reads from, levels of expression, and output location
#
# Output: FASTA files representing reads generated 
#
# Syntax: perl equal_allele_reads.pl single 50 0 ../Dmel/chromFa/chr2L.fa ../Dmel/chromFa/chr2L_alt.fa ../McManus/sim_regs_chr2L.bed ../mel_mel_data/mRNA-Seq/zhr_z30_exons_expression.txt ../Simulations/equal_allele
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
my$expression = $ARGV[6]; # number of transcripts to be generated, sample with replacement
my$out_dir = $ARGV[7]; 

my$refPos = ""; my$refNeg = ""; my$altPos = ""; my$altNeg = "";
our($chr,$pos,$code,$seq,$i,$j,$k,$reflength,$posNeg,$false,$base,$reg,$start,$stop,$begin,$end,$gene,$level,$readRef,$readCons,$sampleRef,$elementsRef,$test);
our(@elements,@parts,@sample);
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

# define expression values for each gene
open(EXP,"$expression") or die "\nError opening $expression\n";
while(<EXP>)    
{
  chomp;
  @elements = split("\t",$_); 
  @parts = split("_",$elements[0]);
  $gene = $parts[0];
  $reg = $parts[1].",".$parts[2];
  $level = $elements[1];
  unless($level == 0){$expr{$gene}{$reg} = $level;}
}
close EXP;

# single end simulation
if($type eq "single")
{
  open(OUT,">$out_dir/chr$chr\_$type\_bp$length\_error$error\_equal_allele.fa") or die "\nError writing to $out_dir/chr$chr\_$type\_bp$length\_error$error\_equal_allele.fa\n";
  foreach $gene (keys %regs)
  {
    #print "\nIn gene loop...";
    foreach $reg (keys %{$regs{$gene}})
    {
      #print "\nGene:$gene\tReg:$reg\tExpression:$expr{$gene}{$reg}\n";
      if($expr{$gene}{$reg} && ($regs{$gene}{$reg}[1]-$regs{$gene}{$reg}[0]) >= $length)
      { 
        #print "\nGene:$gene\tReg:$reg\tStart:$regs{$gene}{$reg}[0]\tStop:$regs{$gene}{$reg}[1]\n";
        #print "Expression:$expr{$gene}{$reg}\n";
        # choose starting points for stochastic read generation

        @elements = ($regs{$gene}{$reg}[0]..($regs{$gene}{$reg}[1]-$length));
        #print "@elements\n";
        $elementsRef = \@elements;
        #$test = int($expr{$gene}{$reg}/4+0.5);
        #print "\nGene: $gene; Region: $reg; Expression value: $expr{$gene}{$reg}; Divided: $test\n";
        #exit;
        $sampleRef = &selection_sample(int($expr{$gene}{$reg}/4+0.5),$elementsRef);
        @sample = @{$sampleRef};
        #print "@sample\n";

        #print "Taken sample...\n";

        # positive orientation
        $k = 1;
        foreach $i (@{$sampleRef})
        {
          if(!$i){die "What the eff is going on?!?!?! Gene:$gene\tReg:$reg\n";}
          #print "i:$i\n";
          #print "The read starting at $i is being generated...\n";

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

          print OUT "\>HWI_$chr\_$gene\_$begin\_$end\_$k\_refPos\n";
          print OUT "$readRef\n";
          print OUT "\>HWI_$chr\_$gene\_$begin\_$end\_$k\_altPos\n";
          print OUT "$readCons\n";
          $k += 1;

        } 
        #exit;

        # choose starting points for stochastic read generation

        @elements = (($reflength-$regs{$gene}{$reg}[1]) .. ($reflength-$regs{$gene}{$reg}[0]-$length));
        $elementsRef = \@elements;
        $sampleRef = &selection_sample(int($expr{$gene}{$reg}/4+0.5),$elementsRef);
        
        # negative orientation
        $k = 1;
        foreach $i (@{$sampleRef})
        {
          #print "i:$i\n";
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
          #print "Begin:$begin\tEnd:$end\n";

          print OUT "\>HWI_$chr\_$gene\_$begin\_$end\_$k\_refNeg\n";
          print OUT "$readRef\n";
          print OUT "\>HWI_$chr\_$gene\_$begin\_$end\_$k\_altNeg\n";
          print OUT "$readCons\n";
          $k += 1;
        }
      }
    }
  }
}
close OUT;

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

sub selection_sample
{
  my($num,$array) = @_;
  my@result;
  my$loop=0;
  my$r;
  while($loop < $num) 
  {
    $r = floor(rand()*scalar(@$array));
    push @result,$array->[$r];
    $loop++;
  }
  return \@result;
}
