#!/usr/bin/perl

###################################################################################
#
# 07/05/2010
#
# analysisPipe.pl
#
# Purpose: analyze gene-specific expression from parents and hybrids to assign 
#          types of cis- and trans-regulatory divergence. The first step will be to
#          filter genes that have less than 20 total coverage. It will then call R scripts 
#          that are responsible for identifying significant differental expression between
#          species. Using this information, the cis- and trans-regulatory divergence
#          assignment will be given to each gene. It will also assess the mode of inheritance 
# 
# Input: 3 files containing the allele count information from mixed parents and hybrids. These will 
#        be combined and exported to R to perform analyses. The data from the analyses will then 
#        be re-imported for cis- and trans-regulatory divergence assignment.
#
# Output: a file containing summary information on each gene, including p-values, CIs, estimates, 
#         raw data, cis-trans reg div assign, mode of inheritance, and whatever else  
#
# E.g. perl analysisPipe.pl <species1> <species2> <geneHyb1> <genePar>
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
  my$s1 = $ARGV[0]; # species 1
  my$s2 = $ARGV[1]; # species 2
  my$hyb = $ARGV[2]; # species 1 X species 2
  my$par = $ARGV[3]; # mixed parents
  my$out_dir = $ARGV[4]; # output directory
  my$ctFilter = $ARGV[5]; # s1 & s2 counts >= ctFilter
  my$pcut = $ARGV[6]; # p-value cutoff (i.e. 0.05)

  my$gene; my$s1ct; my$s2ct; my$tot; # the counts should have already been adjusted in a previous script
  my@elements;
  my%hyb; my%par; my%genes;

  # read in the data from both hybrid crosses and mixed parents

  open(HYB,"$hyb") or die "Can't open $hyb for reading!";
  LINE:while(<HYB>) 
  {
    chomp;
    if($_ =~ /^gene/){next LINE;}
    else
    { 
      @elements = split(/\s/,$_);
      $gene = $elements[0];
      $s1ct = $elements[4];
      $s2ct = $elements[5];
      $tot = $s1ct + $s2ct;
      $genes{$gene} = 1;
      $hyb{$gene} = [$s1ct,$s2ct,$tot];
    }
  } 
  close HYB;

  open(PAR,"$par") or die "Can't open $par for reading!";
  LINE:while(<PAR>) 
  {
    chomp;
    if($_ =~ /^gene/){next LINE;}
    else
    { 
      @elements = split(/\s/,$_);
      $gene = $elements[0];
      $s1ct = $elements[4];
      $s2ct = $elements[5];
      $tot = $s1ct + $s2ct;
      $genes{$gene} = 1;
      $par{$gene} = [$s1ct,$s2ct,$tot];
    }
  } 
  close PAR;

  # Now begin filtering genes with fewer than $ctFilter total counts 

  open(OUT,">$out_dir\/$s1\_$s2.counts.txt") or die "Error writing to $out_dir\/$s1\_$s2.counts.txt\n";
  print OUT "gene\thybS1\thybS2\thybTot\tparS1\tparS2\tparTot";

  foreach $gene (keys %genes)
  {
    if($hyb{$gene} && $par{$gene})
    {
      if($hyb{$gene}[0] >= $ctFilter && $hyb{$gene}[1] >= $ctFilter && $par{$gene}[0] >= $ctFilter && $par{$gene}[1] >= $ctFilter)
      {
        print OUT "\n$gene\t$hyb{$gene}[0]\t$hyb{$gene}[1]\t$hyb{$gene}[2]";
        print OUT "\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]";
      }
    }
  }
  close OUT;
  undef %hyb; undef %par; undef %genes; # clear hashes 

  # Export the data to R for binomial and Fisher's exact tests for cis-trans regulatory divergence

  system "R --no-save < ..\/Rscripts/cis_transMelXSec.R $out_dir\/$s1\_$s2.counts.txt $out_dir\/$s1\_$s2.cis-trans.txt > ./trash";

  # Perform cis-trans regulatory divergence assignment

  my%hyb_par;

  open(IN,"$out_dir\/$s1\_$s2.cis-trans.txt") or die "Can't open $out_dir\/$s1\_$s2.cis-trans.txt for reading!";
  HERE: while(<IN>) 
  {
    chomp;
    if($_ =~ /^gene/){next HERE;}
    else
    {
      @elements = split(/\s/,$_);
      $gene = $elements[0];
      $genes{$gene} = 1;

      # these are the corrected p-values according the FDR method (q-values) for the binomial tests    
      $hyb{$gene} = [$elements[1],$elements[2],$elements[8]]; # [species 1 count, species 2 count, q-val]
      $par{$gene} = [$elements[11],$elements[12],$elements[18]];

      # these are the corrected p-values according to the FDR method (q-values) for the Fisher's exact tests 
      $hyb_par{$gene} = $elements[25];
    }
  } 
  close IN;
  
  sub log2 
  {
    my$n = shift;
    return (log($n)/log(2));
  }
  sub proton
  {
    my$n = shift;
    if($n > 0){return 1;}
    elsif($n < 0){return -1;}
    else{return 0;}
  }

  open(OUT,">$out_dir\/$s1\_$s2.divergence.txt") or die "Error writing to $out_dir\/$s1\_$s2.divergence.txt\n";
  print OUT "gene\thybS1\thybS2\tqvalHyb\tparS1\tparS2\tqvalPar\tqvalHyb_Par\tclass";
  foreach $gene (keys %genes) 
  {
    if($hyb{$gene}[2]<$pcut && $par{$gene}[2]<$pcut && $hyb_par{$gene}<$pcut) # cis+trans & cisXtrans
    {
      if((proton(log2($hyb{$gene}[0]/$hyb{$gene}[1])) == proton(log2($par{$gene}[0]/$par{$gene}[1])))) # cis+trans        
      {
        print OUT "\n$gene\t$hyb{$gene}[0]\t$hyb{$gene}[1]\t$hyb{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb_par{$gene}\tcis+trans";  
      }
      elsif((proton(log2($hyb{$gene}[0]/$hyb{$gene}[1])) != proton(log2($par{$gene}[0]/$par{$gene}[1])))) # cisXtrans
      {
        print OUT "\n$gene\t$hyb{$gene}[0]\t$hyb{$gene}[1]\t$hyb{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb_par{$gene}\tcisXtrans";
      }
    }
    elsif($hyb{$gene}[2]<$pcut && $par{$gene}[2]<$pcut && $hyb_par{$gene}>=$pcut) # cis only
    {
      print OUT "\n$gene\t$hyb{$gene}[0]\t$hyb{$gene}[1]\t$hyb{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb_par{$gene}\tcis";
    }
    elsif($hyb{$gene}[2]>=$pcut && $par{$gene}[2]<$pcut && $hyb_par{$gene}<$pcut) # trans only
    {
      print OUT "\n$gene\t$hyb{$gene}[0]\t$hyb{$gene}[1]\t$hyb{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb_par{$gene}\ttrans";    
    }
    elsif($hyb{$gene}[2]<$pcut && $par{$gene}[2]>=$pcut && $hyb_par{$gene}<$pcut) # compensatory
    {
      print OUT "\n$gene\t$hyb{$gene}[0]\t$hyb{$gene}[1]\t$hyb{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb_par{$gene}\tcompensatory";    
    }
    elsif($hyb{$gene}[2]>=$pcut && $par{$gene}[2]>=$pcut && $hyb_par{$gene}>=$pcut) # conserved
    {
      print OUT "\n$gene\t$hyb{$gene}[0]\t$hyb{$gene}[1]\t$hyb{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb_par{$gene}\tconserved";    
    }
    else # ambiguous
    {
      print OUT "\n$gene\t$hyb{$gene}[0]\t$hyb{$gene}[1]\t$hyb{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb_par{$gene}\tambiguous";    
    }
  }
  close OUT;
}

