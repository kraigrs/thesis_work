#!/usr/bin/perl

######################################################################################
# 
# 03/16/2012
#
# pileup_SNP_ASE.pl
#
# Purpose: read in pileup from SAMtools to look at ASE for SNPs used to generate reads
# 
# Input: list of SNPs, pileup information         
#
# Output: 
#
######################################################################################

use strict;
use warnings;
use Switch;

my$SNPs = $ARGV[0];
my$pileup = $ARGV[1];
my$read_length = $ARGV[2];
my$output = $ARGV[3];

my$chr; my$pos; my$overlaps1; my$overlaps2; my$overlaps3; my$overlaps; my$ref; my$alt; my$base; 
my$ref_allele; my$alt_allele; my$pattern; my$total; my$total1; my$total2; my$total3; my$err;
my$i; my$j; my$k; my$right; my$left; my$SNPs_left; my$SNPs_right; my$steps_right; my$steps_left; my$num;
my@elements; my@bases; my@quals; my@meta; my@positions;
my%SNPhash; my%neighbors;

open(SNPS,"$SNPs") or die "\nError opening $SNPs\n";
while(<SNPS>)
{
  chomp;
  @elements = split("\t",$_);

  # mel_mel
  $chr = $elements[0];
  $pos = $elements[2]; # 1-based position in SNP BED file since 1-based in pileup
  $base = $elements[4]; # if reference genome was used
  #$base = $elements[3]; # if alternative genome was used

  # mel_sim
  #$chr = $elements[0];
  #$pos = $elements[1]; # 1-based position in SNP BED file since 1-based in pileup
  #$base = $elements[3]; # if reference genome was used
  #$base = $elements[2]; # if alternative genome was used


  $SNPhash{$chr}{$pos} = $base;
}
close SNPS;

###### algorithm to compute the number of neighbors within 50 bases

foreach $chr (keys %SNPhash)
{
  foreach $pos (sort {$a <=> $b} keys %{$SNPhash{$chr}})
  {
    push @positions, $pos;
  }
  
  #$num = scalar(@positions);
  #print "$num\n";
  #print "$chr\n"; foreach(@positions){print "$_\n";} #exit;
  
  if(scalar(@positions) == 1){$neighbors{$chr}{$positions[0]} = 0;}
  else
  {
    for($i = 0; $i < @positions; $i++)
    {
      #print "$i\n";
      if($i == 0)
      {
        $steps_right = 0;
        $k = $i + 1;
        $right = $positions[$k] - $positions[$i];

        while($right < $read_length)
        {
          $steps_right += 1;
          if($k == scalar(@positions) - 1){last;}
          $k += 1;
          $right = $positions[$k] - $positions[$i];  
        }

        $neighbors{$chr}{$positions[$i]} = $steps_right;
      }
      elsif($i == scalar(@positions) - 1)
      {
        $steps_left = 0;
        $j = $i - 1;
        $left = $positions[$i] - $positions[$j];

        while($left < $read_length)
        {
          $steps_left += 1;
          if($j == 0){last;}
          $j -= 1;
          $left = $positions[$i] - $positions[$j];  
        }

        $neighbors{$chr}{$positions[$i]} = $steps_left;       
      }
      else
      {
        $steps_right = 0;
        $steps_left = 0;

        $k = $i + 1;
        $j = $i - 1;

        $right = $positions[$k] - $positions[$i];
        $left = $positions[$i] - $positions[$j];

        while($right < $read_length)
        {
          $steps_right += 1;
          if($k == scalar(@positions) - 1){last;}
          $k += 1;
          $right = $positions[$k] - $positions[$i];  
        }

        while($left < $read_length)
        {
          $steps_left += 1;
          if($j == 0){last;}
          $j -= 1;
          $left = $positions[$i] - $positions[$j];  
        }

        $neighbors{$chr}{$positions[$i]} = $steps_right + $steps_left;
      }
    }
  }
  splice @positions;

  #foreach $chr (keys %neighbors)
  #{
  #  foreach $pos (sort {$a <=> $b} keys %{$neighbors{$chr}})
  #  {
  #    print "$chr\t$pos\t$neighbors{$chr}{$pos}\n";
  #  }
  #}
  #exit;
  #print "$chr\n";
}

open(OUT,"> $output") or die "Error writing to $output\n";
print OUT "chr\tpos\ttotal_overlap\tref_allele\talt_allele\terrors\tneighbors";

open(PILE,"$pileup") or die "\nError opening $pileup\n";
while(<PILE>)
{
  chomp;
  @elements = split(/\s+/,$_);
  if($elements[0] =~ /\|/)
  {
    @meta = split(/\|/,$elements[0]);
    $chr = $meta[0];
  }
  else{$chr = $elements[0];}

  $pos = $elements[1];
  $ref_allele = uc($elements[2]);

  $total1 = $elements[3];
  #$total2 = $elements[6];
  #$total3 = $elements[9];
  #$total = $total1+$total2+$total3;
  $total = $total1;

  $overlaps1 = $elements[4];
  #$overlaps2 = $elements[7];
  #$overlaps3 = $elements[10];
  #$overlaps = $overlaps1.$overlaps2.$overlaps3;
  $overlaps = $overlaps1;

  #$total = $elements[3];
  #$overlaps = $elements[4];

  if($SNPhash{$chr}{$pos})
  {
    $ref = 0;
    $alt = 0;
    $err = 0;

    #$alt_allele = &altAllele($ref_allele,$SNPhash{$chr}{$pos});
    $alt_allele = $SNPhash{$chr}{$pos};

    $pattern = $alt_allele.lc($alt_allele);
    #print "\nPattern: [$pattern]\n"; exit;

    @bases = split("",$overlaps);
    foreach(@bases)
    {
      if(/[\.\,]/){$ref += 1;}
      elsif(/[$pattern]/){$alt += 1;}
      elsif(/[ACGTNacgtn]/){$err += 1;}
    }
 
    #if($total >= 20 && $err == 0){print OUT "\n$chr\t$pos\t$total\t$ref\t$alt";}
    print OUT "\n$chr\t$pos\t$total\t$ref\t$alt\t$err\t$neighbors{$chr}{$pos}";
  }
}
close PILE;
close OUT;

sub altAllele
{
  my($n,$amb) = @_;
  my$var;
  switch ($amb)
  {
    case "R"
    {
      if($n eq "A"){$var = "G";}
      elsif($n eq "G"){$var = "A";}
      else{die "Something is seriously wrong here!\n";}
    }
    case "Y"
    {
      if($n eq "C"){$var = "T";}
      elsif($n eq "T"){$var = "C";}
      else{die "Something is seriously wrong here!\n";}
    }
    case "M"
    {
      if($n eq "A"){$var = "C";}
      elsif($n eq "C"){$var = "A";}
      else{die "Something is seriously wrong here!\n";}
    }
    case "K" 
    {
      if($n eq "G"){$var = "T";}
      elsif($n eq "T"){$var = "G";}
      else{die "Something is seriously wrong here!\n";}
    }
    case "S"
    {
      if($n eq "G"){$var = "C";}
      elsif($n eq "C"){$var = "G";}
      else{die "Something is seriously wrong here!\n";}
    }
    case "W"
    {
      if($n eq "A"){$var = "T";}
      elsif($n eq "T"){$var = "A";}
      else{die "Something is seriously wrong here!\n";}
    }
    else{die "Something is seriously wrong here!\n";}
  }
  return $var;
}
