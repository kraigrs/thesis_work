#!/usr/bin/perl

###################################################################################
#
# 07/27/2011
#
# imprinting_pvals.pl
#
# Purpose: generate a list of pvals for each region on eaech chromosome for enrichment
# 
# Input: for each chromosome: 1). list of all genes
#                             2). original imprinting list, in window format
#                             3). estimated null distribution from distance method
#                             4). estimated null distribution from gene cluster method
#
# Output: for each chromosome: 1). pvalue of enrichment of each window from distance method
#                              2). pvalue of enrichment of each window from gene cluster method
#
# Usage: perl imprinting_pvals.pl <const> <imprinting> <null>
#        <permutations>    ==> number of permuted datasets to generate (~500?)
#        <# genes to draw> ==> number of genes to draw from original dataset
#        <const>           ==> list of constitutive genes
#
#  E.g. perl imprinting_pvals.pl ../McManus/constitutive_genes.txt ../mel_mel_data/imprinting ../mel_mel_data/null_dist 500000 10000 50 2L
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
  my$const = $ARGV[0];
  my$imp_dir = $ARGV[1];
  my$null_dir = $ARGV[2];
  my$wind = $ARGV[3];
  my$step = $ARGV[4];
  my$k = $ARGV[5];
  my$chr = $ARGV[6];

  our($start,$stop,$imps,$tot,$coord,$i,$ref,$sum,$perms,$pval,$lengthwind,$lengthelem,$j);
  our(@elements,@windows);
  our(%geneHash,%imprint_dist,%null_dist,%imprint_gene,%null_gene);

  sub num_greater
  {
    my($val,$ref) = @_;
    my@a = @{$ref};
    my$count = 0;
    my$item;
    foreach $item (@a)
    {
      if($item >= $val){$count += 1;}
    }
    return $count;
  }

  #gene hash
  #open(CONST,"$const") or die "Can't open $const for reading!";
  #while(<CONST>) 
  #{
  #  chomp;
  #  if($_ =~ /^#/){next;}
  #  else
  #  { 
  #    @elements = split(/\s+/,$_);
  #    $chr = $elements[1];
  #    if($chr)
  #    {
  #      if($chr ne "-")
  #	 {
  #        $geneHash{$chr} = 1;
  #      }
  #    }
  #  }
  #}
  #close CONST;

  $geneHash{$chr} = 1;

  foreach $chr (keys %geneHash)
  {
    # distance method data
    open(IMP,"$imp_dir\/chr$chr\/zhr_z30\.chr$chr\_wind$wind\_step$step\.imprinting.txt") or die "Can't open $imp_dir\/chr$chr\/zhr_z30\.chr$chr\_wind$wind\_step$step\.imprinting.txt for reading!";
    while(<IMP>) 
    {
      chomp;
      if($_ =~ /left_bound/){next;}
      else
      { 
        @elements = split(/\s+/,$_);
        $start = $elements[0];
        $stop = $elements[1];
        $coord = $start."_".$stop;
        $tot = $elements[2];
        $imps = $elements[3];
        $imprint_dist{$chr}{$coord} = [$imps,$tot];
      }
    }
    close IMP;
    
    $j = 0;
    open(NULL,"$null_dir\/chr$chr\/zhr_z30\.chr$chr\_wind$wind\_step$step\.imprinting.txt") or die "Can't open $null_dir\/chr$chr\/zhr_z30\.chr$chr\_wind$wind\_step$step\.imprinting.txt for reading!";
    while(<NULL>) 
    {
      $j += 1;
      chomp;
      if($_ =~ /^\s+$/){next;}
      elsif($_ =~ /coordinates/)
      {
        @windows = split(/\s+/,$_);
        shift @windows;
      }
      elsif($_ =~ /totals/){next;}
      else
      { 
        @elements = split(/\s+/,$_);
        shift @elements;
        #$lengthwind = scalar(@windows);
        #$lengthelem = scalar(@elements);
        #print "\nfirst element:\{$elements[1]\}\twindow:\{$windows[1]\}\n\n"; 
        #print "length windows:\{$lengthwind\}\tlength elements:\{$lengthelem\}\n\n"; 
        #exit;
        for($i=0;$i<@elements;$i++)
	{ 
          push @{ $null_dist{$chr}{$windows[$i]} },$elements[$i];
        }
      }
    }
    close NULL;
    #$perms = scalar(@{$null_dist{"4"}{"0_500000"}});
    #print "\nj:$j\tperms:$perms\n\n";

    # print pvals from dist method
    open(OUT,">$null_dir\/chr$chr\/zhr_z30\.chr$chr\_wind$wind\_step$step\.imprinting_pvals.txt") or die "Error writing to $null_dir\/chr$chr\/zhr_z30\.chr$chr\_wind$wind\_step$step\.imprinting_pvals.txt\n";
    print OUT "chr\tleft_bound\tright_bound\tpval";
    foreach $coord (keys %{ $null_dist{$chr} })
    {
      if($imprint_dist{$chr}{$coord})
      {
        #print "\nhere!\n\n";
        $ref = \@{ $null_dist{$chr}{$coord} };
        $sum = &num_greater($imprint_dist{$chr}{$coord}[0],$ref);
        $perms = scalar(@{ $null_dist{$chr}{$coord} });
        #print "\nsum: $sum\t# permutations: $perms\n\n";
        $pval = ($sum+1)/($perms+1);
        #if($pval == 0){$pval = 1/$perms;} 
        @elements = split(/_/,$coord);
        print OUT "\n$chr\t$elements[0]\t$elements[1]\t$pval";
        #print "\nchr$chr\tcoord:$coord\tsum:$sum\tperms:$perms\tpval:$pval\n\n";
      }
    }
    close OUT;

    # gene cluster method data
    open(IMP,"$imp_dir\/chr$chr\/zhr_z30\.chr$chr\_k$k.imprinting.txt") or die "Can't open $imp_dir\/chr$chr\/zhr_z30\.chr$chr\_k$k.imprinting.txt for reading!";
    while(<IMP>) 
    {
      chomp;
      if($_ =~ /left_bound/){next;}
      else
      { 
        @elements = split(/\s+/,$_);
        $start = $elements[0];
        $stop = $elements[1];
        $coord = $start."_".$stop;
        $tot = $elements[2];
        $imps = $elements[3];
        $imprint_gene{$chr}{$coord} = [$imps,$tot];
      }
    }
    close IMP;

    open(NULL,"$null_dir\/chr$chr\/zhr_z30\.chr$chr\_k$k.imprinting.txt") or die "Can't open $null_dir\/chr$chr\/zhr_z30\.chr$chr\_k$k.imprinting.txt for reading!";
    while(<NULL>) 
    {
      chomp;
      if($_ =~ /^\s+$/){next;}
      elsif($_ =~ /coordinates/)
      {
        @windows = split(/\s+/,$_);
        shift @windows;
      }
      elsif($_ =~ /totals/){next;}
      else
      { 
        @elements = split(/\s+/,$_);
        shift @elements;
        #$lengthwind = scalar(@windows);
        #$lengthelem = scalar(@elements);
        #print "\nfirst element:\{$elements[1]\}\twindow:\{$windows[1]\}\n\n"; 
        #print "length windows:\{$lengthwind\}\tlength elements:\{$lengthelem\}\n\n"; 
        #exit;
        for($i=0;$i<@elements;$i++)
	{ 
          push @{ $null_gene{$chr}{$windows[$i]} },$elements[$i];
        }
      }
    }
    close NULL;

    # print pvals from gene cluster method
    open(OUT,">$null_dir\/chr$chr\/zhr_z30\.chr$chr\_k$k.imprinting_pvals.txt") or die "Error writing to $null_dir\/chr$chr\/zhr_z30\.chr$chr\_k$k.imprinting_pvals.txt\n";
    print OUT "chr\tleft_bound\tright_bound\tpval";
    foreach $coord (keys %{ $null_gene{$chr} })
    {
      if($imprint_gene{$chr}{$coord})
      {
        $ref = \@{ $null_gene{$chr}{$coord} };
        $sum = &num_greater($imprint_gene{$chr}{$coord}[0],$ref);
        $perms = scalar(@{ $null_gene{$chr}{$coord} });
        #print "\nsum: $sum\t# permutations: $perms\n\n";
        $pval = ($sum+1)/($perms+1); 
        @elements = split(/_/,$coord);
        print OUT "\n$chr\t$elements[0]\t$elements[1]\t$pval";
      }
    }
    close OUT;
  }
}
