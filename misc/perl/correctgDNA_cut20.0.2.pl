#!/usr/bin/perl

###################################################################################
#
# 10/01/2010
#
# correctgDNA_cut20.pl
#
# Purpose: combine summary data from parents and hybrids to obtain common gene set and perform transformation 
# 
# Input: 3 files, 1 summary from mixed parents and 1 summary from each reciprocal cross
#
# Output: a list of common genes between sets that have been adjusted to species-specific sequencing depth and mapping enrichment 
#
# E.g. perl correctgDNA_cut20.pl zhr z30 ../mel_mel_data/mRNA-Seq/zhr/zhr.mRNA.mosaik.geneout.txt 16464075 ../mel_mel_data/mRNA-Seq/z30/z30.mRNA.mosaik.geneout.txt 21806797 ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mRNA.mosaik.geneout.txt 31432754 ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mRNA.mosaik.geneout.txt 31439998 ../mel_mel_data/Resequencing/zhr/zhr.gDNA.mosaik.geneout.txt 42743562 ../mel_mel_data/Resequencing/z30/z30.gDNA.mosaik.geneout.txt 25863911 ../mel_mel_data mosaik.geneout.txt
#
#      <s1>       ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <s2>       ==> indicate name of 2nd species aligned to
#      <parfile> ==> file with parent counts
#      <hyb1file  ==> file with s1 X s2 hybrid counts
#      <hyb1depth>  ==> number of original sequences in hyb 1
#      <hyb2file> ==> file with s2 X s1 hybrid counts
#      <hyb1depth>  ==> number of original sequences in hyb 1
#      <s1gDNA> ==> file with gDNA counts from s1
#      <s1gDNA_depth> ==> gDNA depth for s1
#      <s2gDNA> ==> file with gDNA counts from s2
#      <s2gDNA_depth> ==> gDNA depth for s1
#      <out_dir>   ==> directory to send output
#
###################################################################################

use strict;
use warnings;
use POSIX;

#use Statistics::Descriptive;
#$stat = Statistics::Descriptive::Full->new();
#$stat->add_data(1,2,3,4); $mean = $stat->mean();
#$var  = $stat->variance();
#$tm   = $stat->trimmed_mean(.25);
#$Statistics::Descriptive::Tolerance = 1e-10;

main();
sub main
{
  my$species1 = $ARGV[0]; 
  my$species2 = $ARGV[1];
  my$parfile = $ARGV[2];
  my$hyb1file = $ARGV[3];
  my$hyb1depth = $ARGV[4];
  my$hyb2file = $ARGV[5];
  my$hyb2depth = $ARGV[6];
  my$s1gDNAfile = $ARGV[7];
  my$s1gDNAdepth = $ARGV[8];
  my$s2gDNAfile = $ARGV[9];
  my$s2gDNAdepth = $ARGV[10];
  my$out_dir = $ARGV[11];
  my$suffix = $ARGV[12];

  my$line; my$locus; my$s1ct; my$s2ct; my$both; 
  my$factorHyb1; my$factorHyb2; my$factorS1gDNA; my$factorS2gDNA;
  my$correctS1; my$correctS2;
  my@elements;  
  my%loci; my%par; my%hyb1; my%hyb2; my%s1gDNA; my%s2gDNA;

  # parfile
  open(FILE,"$parfile") or die "Can't open $parfile for reading!";
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
        $par{$locus} = [$s1ct,$s2ct,$both]; 
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

  # s1gDNAfile
  open(FILE,"$s1gDNAfile") or die "Can't open $s1gDNAfile for reading!";
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
      $loci{$locus} = 1;
      $s1gDNA{$locus} = [$s1ct,$s2ct,$both];
    }
  }
  close FILE;

  # s2gDNAfile
  open(FILE,"$s2gDNAfile") or die "Can't open $s2gDNAfile for reading!";
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
      $loci{$locus} = 1;
      $s2gDNA{$locus} = [$s1ct,$s2ct,$both]; 
    }
  }
  close FILE;

  # Determine multiplicative factors for sequencing depth and mapping bias adjustment 

  # sequencing depth adjustment
  
  if($hyb1depth == $hyb2depth){$factorHyb1 = 1; $factorHyb2 = 1;}
  elsif($hyb1depth > $hyb2depth){$factorHyb1 = $hyb2depth/$hyb1depth; $factorHyb2 = 1;}
  elsif($hyb1depth < $hyb2depth){$factorHyb1 = 1; $factorHyb2 = $hyb1depth/$hyb2depth;}

  if($s1gDNAdepth == $s2gDNAdepth){$factorS1gDNA = 1; $factorS2gDNA = 1;}
  elsif($s1gDNAdepth > $s2gDNAdepth){$factorS1gDNA = $s2gDNAdepth/$s1gDNAdepth; $factorS2gDNA = 1;}
  elsif($s1gDNAdepth < $s2gDNAdepth){$factorS1gDNA = 1; $factorS2gDNA = $s1gDNAdepth/$s2gDNAdepth;}

  #print "\nfactorHyb1 = $factorHyb1\nfactorHyb2 = $factorHyb2\n";
  #print "factorS1gDNA = $factorS1gDNA\nfactorS2gDNA = $factorS2gDNA\n";

  open(PARHYB1,">$out_dir/$species1\_$species2.parHyb1_gDNA_cut20.$suffix") or die "Error writing to $out_dir/$species1\_$species2.parHyb1_gDNA_cut20.$suffix\n";
  open(PARHYB2,">$out_dir/$species1\_$species2.parHyb2_gDNA_cut20.$suffix") or die "Error writing to $out_dir/$species1\_$species2.parHyb2_gDNA_cut20.$suffix\n";
  open(HYB1HYB2,">$out_dir/$species1\_$species2.hyb1hyb2_gDNA_cut20.$suffix") or die "Error writing to $out_dir/$species1\_$species2.hyb1hyb2_gDNA_cut20.$suffix\n";

  print PARHYB1  "locus\tPar_$species1\tPar_$species2\tPar_both\tHyb1_$species1\tHyb1_$species2\tHyb1_both\tcorr_$species1\tcorr_$species2";
  print PARHYB2  "locus\tPar_$species1\tPar_$species2\tPar_both\tHyb2_$species1\tHyb2_$species2\tHyb2_both\tcorr_$species1\tcorr_$species2";
  print HYB1HYB2 "locus\tHyb1_$species1\tHyb1_$species2\tHyb1_both\tHyb2_$species1\tHyb2_$species2\tHyb2_both";

  #print "locus\tcorrectS1\tcorrectS2";
  foreach $locus (keys %loci)   
  {
    if($s1gDNA{$locus} && $s2gDNA{$locus} && $s1gDNA{$locus}[0] + $s2gDNA{$locus}[1] >= 20)
    {
      # gene-specific correction
      if($s1gDNA{$locus}[0]*$factorS1gDNA > $s2gDNA{$locus}[1]*$factorS2gDNA)
      {
        $correctS1 = ($s2gDNA{$locus}[1]*$factorS2gDNA)/($s1gDNA{$locus}[0]*$factorS1gDNA);
        $correctS2 = 1;
      }
      else
      {
        $correctS1 = 1;
        $correctS2 = ($s1gDNA{$locus}[0]*$factorS1gDNA)/($s2gDNA{$locus}[1]*$factorS2gDNA);
      }
      if($par{$locus} && $hyb1{$locus})
      {
        printf PARHYB1 ("\n$locus\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
                        floor($par{$locus}[0]*$correctS1),
                        floor($par{$locus}[1]*$correctS2),
                        floor($par{$locus}[2]*$correctS1+$par{$locus}[2]*$correctS2),
                        floor($hyb1{$locus}[0]*$correctS1),
                        floor($hyb1{$locus}[1]*$correctS2),
                        floor($hyb1{$locus}[2]*$correctS1+$hyb1{$locus}[2]*$correctS2));
        print PARHYB1 "\t$correctS1\t$correctS2";
      }
      if($par{$locus} && $hyb2{$locus})
      {
        printf PARHYB2 ("\n$locus\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
                        floor($par{$locus}[0]*$correctS1),
                        floor($par{$locus}[1]*$correctS2),
                        floor($par{$locus}[2]*$correctS1+$par{$locus}[2]*$correctS2),
                        floor($hyb2{$locus}[0]*$correctS1),
                        floor($hyb2{$locus}[1]*$correctS2),
                        floor($hyb2{$locus}[2]*$correctS1+$hyb2{$locus}[2]*$correctS2));
        print PARHYB2 "\t$correctS1\t$correctS2";
      }
      if($hyb1{$locus} && $hyb2{$locus})
      {
        printf HYB1HYB2 ("\n$locus\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
                        floor($hyb1{$locus}[0]*$factorHyb1*$correctS1),
                        floor($hyb1{$locus}[1]*$factorHyb1*$correctS2),
                        floor($hyb1{$locus}[2]*$factorHyb1*$correctS1+$hyb1{$locus}[2]*$factorHyb1*$correctS2),
                        floor($hyb2{$locus}[0]*$factorHyb2*$correctS1),
                        floor($hyb2{$locus}[1]*$factorHyb2*$correctS2),
                        floor($hyb2{$locus}[2]*$factorHyb2*$correctS1+$hyb2{$locus}[2]*$factorHyb2*$correctS2));
      }
    }
  }
  close PARHYB1; close PARHYB2; close HYB1HYB2;
}
