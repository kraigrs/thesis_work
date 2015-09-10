#!/usr/bin/perl

######################################################################################
# 
# 08/21/2012
#
# combine_mpileups_SNP_ASE.pl
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

my$ref = $ARGV[0];
my$alt = $ARGV[1];
my$SNPs = $ARGV[2];
my$ref_pileup = $ARGV[3];
my$alt_pileup = $ARGV[4];
my$read_length = $ARGV[5];
my$output = $ARGV[6];

our($chr,$pos,$overlaps,$base1,$base2,$ref_allele,$alt_allele,$pattern,$total,$i,$j,$k,$right,$left,$SNPs_left,$SNPs_right,$steps_right,$steps_left,$num);
our(@elements,@bases,@positions);
our(%SNPhash,%refHash,%altHash,%neighbors);

open(SNPS,"$SNPs") or die "\nError opening $SNPs\n";
while(<SNPS>)
{
  chomp;
  @elements = split(/\s+/,$_);
  $chr = $elements[0];
  $pos = $elements[1]; # 1-based position in SNP BED file since 1-based in pileup
  $base1 = $elements[2];
  $base2 = $elements[3];

  $SNPhash{$chr}{$pos} = [$base1,$base2];

  $refHash{$chr}{$pos}{'ref'} = 0;
  $refHash{$chr}{$pos}{'alt'} = 0;

  $altHash{$chr}{$pos}{'ref'} = 0;
  $altHash{$chr}{$pos}{'alt'} = 0;

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

print "\nhere\n";

open(PILE1,"$ref_pileup") or die "\nError opening $ref_pileup\n";
while(<PILE1>)
{
  chomp;

  @elements = split(/\s+/,$_);
  $chr = $elements[0];
  $pos = $elements[1];
  $ref_allele = uc($elements[2]);
  $overlaps = $elements[4].$elements[7].$elements[10];

  #print "Chr: $chr\tPos: $pos\n"; exit; 
  
  if($SNPhash{$chr}{$pos})
  {
    $alt_allele = $SNPhash{$chr}{$pos}[1];

    $pattern = $alt_allele.lc($alt_allele);

    @bases = split("",$overlaps);
    foreach(@bases)
    {
      if(/[\.,]/){$refHash{$chr}{$pos}{'ref'} += 1;}
      elsif(/[$pattern]/){$refHash{$chr}{$pos}{'alt'} += 1;}
    }
  }
}
close PILE1;

print "\nhere\n";

open(PILE2,"$alt_pileup") or die "\nError opening $alt_pileup\n";
while(<PILE2>)
{
  chomp;

  @elements = split(/\s+/,$_);
  $chr = $elements[0];
  $pos = $elements[1];
  $alt_allele = uc($elements[2]);
  $overlaps = $elements[4].$elements[7].$elements[10];
  
  if($SNPhash{$chr}{$pos})
  {
    $ref_allele = $SNPhash{$chr}{$pos}[0];

    $pattern = $ref_allele.lc($ref_allele);

    @bases = split("",$overlaps);
    foreach(@bases)
    {
      if(/[\.,]/){$altHash{$chr}{$pos}{'alt'} += 1;}
      elsif(/[$pattern]/){$altHash{$chr}{$pos}{'ref'} += 1;}
    }
  }
}
close PILE2;

print "\nhere\n";

open(OUT,"> $output") or die "Error writing to $output\n";
print OUT "chr\tpos\t$ref\_ref_allele\t$ref\_alt_allele\t$alt\_ref_allele\t$alt\_alt_allele\tneighbors\n";

foreach $chr (keys %SNPhash)
{
  foreach $pos (sort {$a <=> $b} keys %{$SNPhash{$chr}})
  {

    #print "$chr\t$pos\t$refHash{$chr}{$pos}{ref}\t$refHash{$chr}{$pos}{alt}\t$altHash{$chr}{$pos}{ref}\t$altHash{$chr}{$pos}{alt}\t$neighbors{$chr}{$pos}\n";

    #if(!$refHash{$chr}{$pos}{ref}||!$refHash{$chr}{$pos}{alt}||!$altHash{$chr}{$pos}{ref}||!$altHash{$chr}{$pos}{alt}||!$neighbors{$chr}{$pos})
    #{
    #  print "$chr\t$pos\n";
    #}

    #else{print "$chr\t$pos\t$refHash{$chr}{$pos}{ref}\t$refHash{$chr}{$pos}{alt}\t$altHash{$chr}{$pos}{ref}\t$altHash{$chr}{$pos}{alt}\n";}

    print OUT "$chr\t$pos\t$refHash{$chr}{$pos}{ref}\t$refHash{$chr}{$pos}{alt}\t$altHash{$chr}{$pos}{ref}\t$altHash{$chr}{$pos}{alt}\t$neighbors{$chr}{$pos}\n";
  }
}

close OUT;
