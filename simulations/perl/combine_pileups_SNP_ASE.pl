#!/usr/bin/perl

######################################################################################
# 
# 04/13/2012
#
# combine_pileups_SNP_ASE.pl
#
# Purpose: read in pileups generated from aligning to separate genomes and combine ASE measures
# 
# Input: list of SNPs, pileup information         
#
# Output: 
#
######################################################################################

use strict;
use warnings;
use Switch;

my$ref = $ARGV[0];
my$alt = $ARGV[1];
my$SNPs = $ARGV[2];
my$ref_pileup = $ARGV[3];
my$alt_pileup = $ARGV[4];
my$output = $ARGV[5];

our($chr,$pos,$overlaps,$base,$ref_allele,$alt_allele,$pattern,$total,$i,);
our(@elements,@bases,@quals,@positions);
our(%SNPhash,%refHash,%altHash,%genome,%neighbors);

my$found = 0; my$seq = "";

#print "\nBuilding reference hash...\n";

open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  #print "Here!\n";
  if(/^>(chr\S+)$/ && $found == 1)
  {
    $genome{$chr} = $seq;

    $seq = "";
    $chr = $1;
  }
  elsif(/^>(chr\S+)$/)
  {
    $seq = "";
    $chr = $1;
    $found = 1;
  }
  elsif($found == 1)
  {
    $seq = $seq.uc($_);
  }
}
$genome{$chr} = $seq;
close REF;

open(SNPS,"$SNPs") or die "\nError opening $SNPs\n";
while(<SNPS>)
{
  chomp;
  @elements = split("\t",$_);
  $chr = $elements[0];
  $pos = $elements[2]; # 1-based position in SNP BED file since 1-based in pileup
  $base = $elements[3];
  $SNPhash{$chr}{$pos} = $base;

  $refHash{$chr}{$pos}{'ref'} = 0;
  $refHash{$chr}{$pos}{'alt'} = 0;

  $altHash{$chr}{$pos}{'ref'} = 0;
  $altHash{$chr}{$pos}{'alt'} = 0;

}
close SNPS;

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

open(PILE1,"$ref_pileup") or die "\nError opening $ref_pileup\n";
while(<PILE1>)
{
  chomp;

  @elements = split("\t",$_);
  $chr = $elements[0];
  $pos = $elements[1];
  $ref_allele = uc($elements[2]);
  $total = $elements[3];
  $overlaps = $elements[4];

  #print "Chr: $chr\tPos: $pos\n"; exit; 
  
  if($SNPhash{$chr}{$pos})
  {
    #print "Here\n"; exit;
    #$alt_allele = &altAllele($ref_allele,$SNPhash{$chr}{$pos});
    $alt_allele = $SNPhash{$chr}{$pos};

    $pattern = $alt_allele.lc($alt_allele);
    #print "\nPattern: [$pattern]\n"; exit;

    @bases = split("",$overlaps);
    foreach(@bases)
    {
      if(/[\.,]/){$refHash{$chr}{$pos}{'ref'} += 1;}
      elsif(/[$pattern]/){$refHash{$chr}{$pos}{'alt'} += 1;}
    }
  }
}
close PILE1;

open(PILE2,"$alt_pileup") or die "\nError opening $alt_pileup\n";
while(<PILE2>)
{
  chomp;

  @elements = split("\t",$_);
  $chr = $elements[0];
  $pos = $elements[1];
  $ref_allele = uc($elements[2]);
  $total = $elements[3];
  $overlaps = $elements[4];
  
  if($SNPhash{$chr}{$pos})
  {
    #$alt_allele = &altAllele($ref_allele,$SNPhash{$chr}{$pos});
    $alt_allele = substr $genome{$chr},$pos-1,1;

    $pattern = $alt_allele.lc($alt_allele);
    #print "\nPattern: [$pattern]\n"; exit;

    @bases = split("",$overlaps);
    foreach(@bases)
    {
      if(/[\.,]/){$altHash{$chr}{$pos}{'alt'} += 1;}
      elsif(/[$pattern]/){$altHash{$chr}{$pos}{'ref'} += 1;}
    }
  }
}
close PILE2;

open(OUT,"> $output") or die "Error writing to $output\n";
print OUT "chr\tpos\tdm3_ref_ref_allele\tdm3_ref_alt_allele\tdm3_alt_ref_allele\tdm3_alt_alt_allele\n";

foreach $chr (keys %SNPhash)
{
  foreach $pos (keys %{$SNPhash{$chr}})
  {
    #print "\n$chr\t$pos\t$refHash{$chr}{$pos}{ref}\t$refHash{$chr}{$pos}{alt}\t$altHash{$chr}{$pos}{ref}\t$altHash{$chr}{$pos}{alt}\n"; exit;
    print OUT "$chr\t$pos\t$refHash{$chr}{$pos}{ref}\t$refHash{$chr}{$pos}{alt}\t$altHash{$chr}{$pos}{ref}\t$altHash{$chr}{$pos}{alt}\n";
  }
}

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
