#!/usr/bin/perl

###################################################################################
#
# 10/01/2010
#
# correctMapBias_quant.pl
#
# Purpose: combine summary data from parents and hybrids to obtain common gene set and perform transformation 
# 
# Input: 4 files, 1 summary from each parent and 1 summary from each reciprocal cross
#
# Output: a list of common genes between sets that have been adjusted to species-specific sequencing depth and mapping enrichment 
#
# E.g. perl processData.pl <species1> <species2> <species1file> <s1depth> <species2file> <s2depth> <hyb1file> <hyb1depth> <hyb2file> <hyb2depth> <out_dir>
#
#      <s1>       ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <s2>       ==> indicate name of 2nd species aligned to
#      <par1file> ==> file with s1 parent counts
#      <par2file> ==> file with s2 parent counts
#      <hyb1file  ==> file with s1 X s2 hybrid counts
#      <hyb2file> ==> file with s2 X s1 hybrid counts
#      <out_dir>   ==> directory to send output
#
###################################################################################

use strict;
use warnings;
use POSIX;

main();
sub main
{
  my$species1 = $ARGV[0]; 
  my$species2 = $ARGV[1];
  my$par1file = $ARGV[2];
  my$par2file = $ARGV[3];
  my$hyb1file = $ARGV[4];
  my$hyb2file = $ARGV[5];
  my$out_dir = $ARGV[6];
  my$suffix = $ARGV[7];

  my$line; my$locus; my$s1ct; my$s2ct; my$both; 
  my$parHyb1s1=0; my$parHyb1s2=0; my$parHyb2s1=0; my$parHyb2s2=0; my$hyb1tot=0; my$hyb2tot=0;
  my$factorParHyb1s1; my$factorParHyb1s2; 
  my$factorParHyb2s1; my$factorParHyb2s2;
  my$factorHyb1; my$factorHyb2;
  my$LBparHyb1s1; my$UBparHyb1s1; my$LBparHyb1s2; my$UBparHyb1s2;
  my$LBparHyb2s1; my$UBparHyb2s1; my$LBparHyb2s2; my$UBparHyb2s2;
  my$parHyb1s1ctsRef; my$parHyb1s2ctsRef;
  my$parHyb2s1ctsRef; my$parHyb2s2ctsRef;
  my@elements; my@parHyb1s1cts; my@parHyb1s2cts; my@parHyb2s1cts; my@parHyb2s2cts; 
  my%loci; my%par1; my%par2; my%hyb1; my%hyb2;

  # Quantization
  # Usage: quantile7(<arrayRef>,<p>) = interpolated quantile
  sub quantile7
  {
    my($arrayRef,$p) = @_;
    my@sorted = sort {$a <=> $b} @{$arrayRef};
    my$N = @sorted;
    my$h = ($N-1)*$p + 1;
    my$Qp = $sorted[floor($h)-1] + ($h - floor($h))*($sorted[floor($h)] - $sorted[floor($h)-1]);
    return $Qp;
  }

  # par1file
  open(FILE,"$par1file") or die "Can't open $par1file for reading!";
  LINE:while(<FILE>) 
  {
    chomp;
    $line = $_;  
    if($line =~ /^(gene|exon).+/){next LINE;}
    elsif($line =~ /^CG.+/)
    { 
      @elements = split(/\s/,$line);
      $locus = $elements[0];
      $s1ct = $elements[1];
      $s2ct = $elements[2];
      $both = $elements[3];
      if($s1ct + $s2ct > 20)
      {
        $loci{$locus} = 1;
        $par1{$locus} = [$s1ct,$s2ct,$both];
      } 
    }
  }
  close FILE;

  # par2file
  open(FILE,"$par2file") or die "Can't open $par2file for reading!";
  LINE:while(<FILE>) 
  {
    chomp;
    $line = $_;  
    if($line =~ /^(gene|exon).+/){next LINE;}
    elsif($line =~ /^CG.+/)
    { 
      @elements = split(/\s/,$line);
      $locus = $elements[0];
      $s1ct = $elements[1];
      $s2ct = $elements[2];
      $both = $elements[3];
      if($s1ct + $s2ct > 20)
      {
        $loci{$locus} = 1;
        $par2{$locus} = [$s1ct,$s2ct,$both]; 
      }
    }
  }
  close FILE;

  # hyb1file
  open(FILE,"$hyb1file") or die "Can't open $hyb1file for reading!";
  LINE:while(<FILE>) 
  {
    chomp;
    $line = $_;  
    if($line =~ /^(gene|exon).+/){next LINE;}
    elsif($line =~ /^CG.+/)
    { 
      @elements = split(/\s/,$line);
      $locus = $elements[0];
      $s1ct = $elements[1];
      $s2ct = $elements[2];
      $both = $elements[3];
      if($s1ct + $s2ct > 20)
      {
        $loci{$locus} = 1;
        $hyb1{$locus} = [$s1ct,$s2ct,$both];
      }
    }
  }
  close FILE;

  # hyb2file
  open(FILE,"$hyb2file") or die "Can't open $hyb2file for reading!";
  LINE:while(<FILE>) 
  {
    chomp;
    $line = $_;  
    if($line =~ /^(gene|exon).+/){next LINE;}
    elsif($line =~ /^CG.+/)
    { 
      @elements = split(/\s/,$line);
      $locus = $elements[0];
      $s1ct = $elements[1];
      $s2ct = $elements[2];
      $both = $elements[3];
      if($s1ct + $s2ct > 20)
      {
        $loci{$locus} = 1;
        $hyb2{$locus} = [$s1ct,$s2ct,$both]; 
      }
    }
  }
  close FILE;

  # Determine accurate counts for each case using common gene sets

  foreach $locus (keys %loci)   
  {
    if($par1{$locus} && $par2{$locus} && $hyb1{$locus})
    {
      push @parHyb1s1cts,$par1{$locus}[0];
      push @parHyb1s2cts,$par2{$locus}[1];
    }
    if($par1{$locus} && $par2{$locus} && $hyb2{$locus})
    {
      push @parHyb2s1cts,$par1{$locus}[0];
      push @parHyb2s2cts,$par2{$locus}[1];
    }
    if($hyb1{$locus} && $hyb2{$locus})
    {
      $hyb1tot += $hyb1{$locus}[0] + $hyb1{$locus}[1];
      $hyb2tot += $hyb2{$locus}[0] + $hyb2{$locus}[1];
    }
  }

  $parHyb1s1ctsRef = \@parHyb1s1cts; $parHyb1s2ctsRef = \@parHyb1s2cts;
  $parHyb2s1ctsRef = \@parHyb2s1cts; $parHyb2s2ctsRef = \@parHyb2s2cts;

  $LBparHyb1s1 = quantile7($parHyb1s1ctsRef,0.025); $UBparHyb1s1 = quantile7($parHyb1s1ctsRef,0.975);
  $LBparHyb1s2 = quantile7($parHyb1s2ctsRef,0.025); $UBparHyb1s2 = quantile7($parHyb1s2ctsRef,0.975);
  $LBparHyb2s1 = quantile7($parHyb2s1ctsRef,0.025); $UBparHyb2s1 = quantile7($parHyb2s1ctsRef,0.975);
  $LBparHyb2s2 = quantile7($parHyb2s2ctsRef,0.025); $UBparHyb2s2 = quantile7($parHyb2s2ctsRef,0.975);

  print "LBparHyb1s1 = $LBparHyb1s1\nUBparHyb1s1 = $UBparHyb1s1\nLBparHyb1s2 = $LBparHyb1s2\nUBparHyb1s2 = $UBparHyb1s2\n\n";
  print "LBparHyb2s1 = $LBparHyb2s1\nUBparHyb2s1 = $UBparHyb2s1\nLBparHyb2s2 = $LBparHyb2s2\nUBparHyb2s2 = $UBparHyb2s2\n\n";

  foreach $locus (keys %loci)   
  {
    if($par1{$locus} && $par2{$locus} && $hyb1{$locus} && 
       $par1{$locus}[0] > $LBparHyb1s1 && $par1{$locus}[0] < $UBparHyb1s1 &&
       $par2{$locus}[1] > $LBparHyb1s2 && $par2{$locus}[1] < $UBparHyb1s2)
    {
      $parHyb1s1 += $par1{$locus}[0]; $parHyb1s2 += $par2{$locus}[1];
    }
    if($par1{$locus} && $par2{$locus} && $hyb2{$locus} && 
       $par1{$locus}[0] > $LBparHyb2s1 && $par1{$locus}[0] < $UBparHyb2s1 &&
       $par2{$locus}[1] > $LBparHyb2s2 && $par2{$locus}[1] < $UBparHyb2s2)
    {
      $parHyb2s1 += $par1{$locus}[0]; $parHyb2s2 += $par2{$locus}[1];
    }
  }

  # Determine multiplicative factors for sequencing depth and mapping bias adjustment 
  
  # sequencing depth and mapping bias adjustment 
  if($parHyb1s1 == $parHyb1s2){$factorParHyb1s1 = 1; $factorParHyb1s2 = 1;}
  elsif($parHyb1s1 > $parHyb1s2){$factorParHyb1s1 = $parHyb1s2/$parHyb1s1; $factorParHyb1s2 = 1;}
  elsif($parHyb1s1 < $parHyb1s2){$factorParHyb1s1 = 1; $factorParHyb1s2 = $parHyb1s1/$parHyb1s2;}

  if($parHyb2s1 == $parHyb2s2){$factorParHyb2s1 = 1; $factorParHyb2s2 = 1;}
  elsif($parHyb2s1 > $parHyb2s2){$factorParHyb2s1 = $parHyb2s2/$parHyb2s1; $factorParHyb2s2 = 1;}
  elsif($parHyb2s1 < $parHyb2s2){$factorParHyb2s1 = 1; $factorParHyb2s2 = $parHyb2s1/$parHyb2s2;}
 
  if($hyb1tot == $hyb2tot){$factorHyb1 = 1; $factorHyb2 = 1;}
  elsif($hyb1tot > $hyb2tot){$factorHyb1 = $hyb2tot/$hyb1tot; $factorHyb2 = 1;}
  elsif($hyb1tot < $hyb2tot){$factorHyb1 = 1; $factorHyb2 = $hyb1tot/$hyb2tot;}

  print "\nfactorParHyb1s1 = $factorParHyb1s1\nfactorParHyb1s2 = $factorParHyb1s2\n";
  print "factorParHyb2s1 = $factorParHyb2s1\nfactorParHyb2s2 = $factorParHyb2s2\n";
  print "factorHyb1 = $factorHyb1\nfactorHyb2 = $factorHyb2\n\n";

  open(PARHYB1,">$out_dir/$species1\_$species2.parHyb1_mapBias_quant.$suffix") or die "Error writing to $out_dir/$species1\_$species2.parHyb1.$suffix\n";
  open(PARHYB2,">$out_dir/$species1\_$species2.parHyb2_mapBias_quant.$suffix") or die "Error writing to $out_dir/$species1\_$species2.parHyb2.$suffix\n";
  open(HYB1HYB2,">$out_dir/$species1\_$species2.hyb1hyb2_mapBias_quant.$suffix") or die "Error writing to $out_dir/$species1\_$species2.hyb1hyb2.$suffix\n";

  print PARHYB1  "locus\tPar_$species1\tPar_$species2\tPar_both\tHyb1_$species1\tHyb1_$species2\tHyb1_both";
  print PARHYB2  "locus\tPar_$species1\tPar_$species2\tPar_both\tHyb2_$species1\tHyb2_$species2\tHyb2_both";
  print HYB1HYB2 "locus\tHyb1_$species1\tHyb1_$species2\tHyb1_both\tHyb2_$species1\tHyb2_$species2\tHyb2_both";

  foreach $locus (keys %loci)   
  {
    if($par1{$locus} && $par2{$locus} && $hyb1{$locus} && 
       $par1{$locus}[0] > $LBparHyb1s1 && $par1{$locus}[0] < $UBparHyb1s1 &&
       $par2{$locus}[1] > $LBparHyb1s2 && $par2{$locus}[1] < $UBparHyb1s2)
    {
      printf PARHYB1 ("\n$locus\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f",
                      $par1{$locus}[0]*$factorParHyb1s1,$par2{$locus}[1]*$factorParHyb1s2,
                      $par1{$locus}[2]*$factorParHyb1s1+$par2{$locus}[2]*$factorParHyb1s2,
                      $hyb1{$locus}[0],$hyb1{$locus}[1],$hyb1{$locus}[2]);
    }
    if($par1{$locus} && $par2{$locus} && $hyb2{$locus} && 
       $par1{$locus}[0] > $LBparHyb2s1 && $par1{$locus}[0] < $UBparHyb2s1 &&
       $par2{$locus}[1] > $LBparHyb2s2 && $par2{$locus}[1] < $UBparHyb2s2)
    {
      printf PARHYB2 ("\n$locus\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f",
                      $par1{$locus}[0]*$factorParHyb2s1,$par2{$locus}[1]*$factorParHyb2s2,
                      $par1{$locus}[2]*$factorParHyb2s1+$par2{$locus}[2]*$factorParHyb2s2,
                      $hyb2{$locus}[0],$hyb2{$locus}[1],$hyb2{$locus}[2]);
    }
    if($hyb1{$locus} && $hyb2{$locus})
    {
      printf HYB1HYB2 ("\n$locus\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f",
                      $hyb1{$locus}[0]*$factorHyb1,$hyb1{$locus}[1]*$factorHyb1,$hyb1{$locus}[2]*$factorHyb1,
                      $hyb2{$locus}[0]*$factorHyb2,$hyb2{$locus}[1]*$factorHyb2,$hyb2{$locus}[2]*$factorHyb2);
    }
  }
  close PARHYB1; close PARHYB2; close HYB1HYB2;
}
