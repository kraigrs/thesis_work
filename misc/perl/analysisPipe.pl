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
# E.g. perl analysisPipe.pl <species1> <species2> <geneHyb1> <geneHyb2> <genePar>
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
  my$s1 = $ARGV[0]; # species 1
  my$s2 = $ARGV[1]; # species 2
  my$hyb1 = $ARGV[2]; # species 1 X species 2
  my$hyb2 = $ARGV[3]; # species 2 X species 1
  my$par = $ARGV[4]; # mixed parents
  my$out_dir = $ARGV[5]; # output directory
  my$ctFilter = $ARGV[6]; # s1 & s2 counts >= ctFilter
  my$pcut = $ARGV[7]; # p-value cutoff (i.e. 0.05)

  my$gene; my$s1ct; my$s2ct; my$tot; # the counts should have already been adjusted in a previous script
  my@elements;
  my%hyb1; my%hyb2; my%par; my%genes;

  # read in the data from both hybrid crosses and mixed parents

  open(HYB1,"$hyb1") or die "Can't open $hyb1 for reading!";
  LINE:while(<HYB1>) 
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
      $hyb1{$gene} = [$s1ct,$s2ct,$tot];
    }
  } 
  close HYB1;

  open(HYB2,"$hyb2") or die "Can't open $hyb2 for reading!";
  LINE:while(<HYB2>) 
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
      $hyb2{$gene} = [$s1ct,$s2ct,$tot];
    }
  } 
  close HYB2;

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
    if($hyb1{$gene} && $par{$gene})
    {
      if($hyb1{$gene}[0] >= $ctFilter && $hyb1{$gene}[1] >= $ctFilter && $par{$gene}[0] >= $ctFilter && $par{$gene}[1] >= $ctFilter)
      {
        print OUT "\n$gene\t$hyb1{$gene}[0]\t$hyb1{$gene}[1]\t$hyb1{$gene}[2]";
        print OUT "\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]";
      }
    }
  }
  close OUT;

  open(OUT,">$out_dir\/$s2\_$s1.counts.txt") or die "Error writing to $out_dir\/$s2\_$s1.counts.txt\n";
  print OUT "gene\thybS1\thybS2\thybTot\tparS1\tparS2\tparTot";

  foreach $gene (keys %genes)
  {
    if($hyb2{$gene} && $par{$gene})
    {
      if($hyb2{$gene}[0] >= $ctFilter && $hyb2{$gene}[1] >= $ctFilter && $par{$gene}[0] >= $ctFilter && $par{$gene}[1] >= $ctFilter)
      {
        print OUT "\n$gene\t$hyb2{$gene}[0]\t$hyb2{$gene}[1]\t$hyb2{$gene}[2]";
        print OUT "\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]";
      }
    }
  }
  close OUT;

  undef %hyb1; undef %hyb2; undef %par; undef %genes; # clear hashes 

  # Export the data to R for binomial and Fisher's exact tests for cis-trans regulatory divergence

  system "R --no-save < ../Rscripts/cis_transMelXSec.R $out_dir/$s1\_$s2.counts.txt $out_dir/$s1\_$s2.cis-trans.txt > ./trash";
  system "R --no-save < ../Rscripts/cis_transMelXSec.R $out_dir/$s2\_$s1.counts.txt $out_dir/$s2\_$s1.cis-trans.txt > ./trash";

  # Perform cis-trans regulatory divergence assignment

  my%hyb1_par; my%hyb2_par; my%hyb1_hyb2; my%genes1; my%genes2;

  open(IN,"$out_dir/$s1\_$s2.cis-trans.txt") or die "Can't open $out_dir/$s1\_$s2.cis-trans.txt for reading!";
  HERE: while(<IN>) 
  {
    chomp;
    if($_ =~ /^gene/){next HERE;}
    else
    {    
      @elements = split(/\s/,$_);
      $gene = $elements[0];
      $genes1{$gene} = 1;

      # these are the corrected p-values according the FDR method (q-values) for the binomial tests    
      $hyb1{$gene} = [$elements[1],$elements[2],$elements[8]]; # [species 1 count, species 2 count, q-val]
      $par{$gene} = [$elements[11],$elements[12],$elements[18]];

      # these are the corrected p-values according to the FDR method (q-values) for the Fisher's exact tests 
      $hyb1_par{$gene} = $elements[25];
    }
  } 
  close IN;

  open(IN,"$out_dir/$s2\_$s1.cis-trans.txt") or die "Can't open $out_dir/$s2\_$s1.cis-trans.txt for reading!";
  HERE: while(<IN>) 
  {
    chomp;
    if($_ =~ /^gene/){next HERE;}
    else
    {    
      @elements = split(/\s/,$_);
      $gene = $elements[0];
      $genes2{$gene} = 1;

      # these are the corrected p-values according the FDR method (q-values) for the binomial tests    
      $hyb2{$gene} = [$elements[1],$elements[2],$elements[8]];
      $par{$gene} = [$elements[11],$elements[12],$elements[18]];

      # these are the corrected p-values according to the FDR method (q-values) for the Fisher's exact tests 
      $hyb2_par{$gene} = $elements[25];
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

  # hyb1 cross
  open(OUT,">$out_dir/$s1\_$s2.divergence.txt") or die "Error writing to $out_dir/$s1\_$s2.divergence.txt\n";
  print OUT "gene\thybS1\thybS2\tqvalHyb\tparS1\tparS2\tqvalPar\tqvalHyb_Par\tclass";
  foreach $gene (keys %genes1) 
  {
    if($hyb1{$gene}[2]<$pcut && $par{$gene}[2]<$pcut && $hyb1_par{$gene}<$pcut) # cis+trans & cisXtrans
    {
      if((proton(log2($hyb1{$gene}[0]/$hyb1{$gene}[1])) == proton(log2($par{$gene}[0]/$par{$gene}[1])))) # cis+trans        
      {
        print OUT "\n$gene\t$hyb1{$gene}[0]\t$hyb1{$gene}[1]\t$hyb1{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb1_par{$gene}\tcis+trans";  
      }
      elsif((proton(log2($hyb1{$gene}[0]/$hyb1{$gene}[1])) != proton(log2($par{$gene}[0]/$par{$gene}[1])))) # cisXtrans
      {
        print OUT "\n$gene\t$hyb1{$gene}[0]\t$hyb1{$gene}[1]\t$hyb1{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb1_par{$gene}\tcisXtrans";
      }
    }
    elsif($hyb1{$gene}[2]<$pcut && $par{$gene}[2]<$pcut && $hyb1_par{$gene}>=$pcut) # cis only
    {
      print OUT "\n$gene\t$hyb1{$gene}[0]\t$hyb1{$gene}[1]\t$hyb1{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb1_par{$gene}\tcis";
    }
    elsif($hyb1{$gene}[2]>=$pcut && $par{$gene}[2]<$pcut && $hyb1_par{$gene}<$pcut) # trans only
    {
      print OUT "\n$gene\t$hyb1{$gene}[0]\t$hyb1{$gene}[1]\t$hyb1{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb1_par{$gene}\ttrans";    
    }
    elsif($hyb1{$gene}[2]<$pcut && $par{$gene}[2]>=$pcut && $hyb1_par{$gene}<$pcut) # compensatory
    {
      print OUT "\n$gene\t$hyb1{$gene}[0]\t$hyb1{$gene}[1]\t$hyb1{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb1_par{$gene}\tcompensatory";    
    }
    elsif($hyb1{$gene}[2]>=$pcut && $par{$gene}[2]>=$pcut && $hyb1_par{$gene}>=$pcut) # conserved
    {
      print OUT "\n$gene\t$hyb1{$gene}[0]\t$hyb1{$gene}[1]\t$hyb1{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb1_par{$gene}\tconserved";    
    }
    else # ambiguous
    {
      print OUT "\n$gene\t$hyb1{$gene}[0]\t$hyb1{$gene}[1]\t$hyb1{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb1_par{$gene}\tambiguous";    
    }
  }
  close OUT;

  # hyb2 cross
  open(OUT,">$out_dir/$s2\_$s1.divergence.txt") or die "Error writing to $out_dir/$s2\_$s1.divergence.txt\n";
  print OUT "gene\thybS1\thybS2\tqvalHyb\tparS1\tparS2\tqvalPar\tqvalHyb_Par\tclass";
  foreach $gene (keys %genes2) 
  {
    if($hyb2{$gene}[2]<$pcut && $par{$gene}[2]<$pcut && $hyb2_par{$gene}<$pcut) # cis+trans & cisXtrans
    {
      if((proton(log2($hyb2{$gene}[0]/$hyb2{$gene}[1])) == proton(log2($par{$gene}[0]/$par{$gene}[1])))) # cis+trans        
      {
        print OUT "\n$gene\t$hyb2{$gene}[0]\t$hyb2{$gene}[1]\t$hyb2{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb2_par{$gene}\tcis+trans";  
      }
      elsif((proton(log2($hyb2{$gene}[0]/$hyb2{$gene}[1])) != proton(log2($par{$gene}[0]/$par{$gene}[1])))) # cisXtrans
      {
        print OUT "\n$gene\t$hyb2{$gene}[0]\t$hyb2{$gene}[1]\t$hyb2{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb2_par{$gene}\tcisXtrans";
      }
    }
    elsif($hyb2{$gene}[2]<$pcut && $par{$gene}[2]<$pcut && $hyb2_par{$gene}>=$pcut) # cis only
    {
      print OUT "\n$gene\t$hyb2{$gene}[0]\t$hyb2{$gene}[1]\t$hyb2{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb2_par{$gene}\tcis";
    }
    elsif($hyb2{$gene}[2]>=$pcut && $par{$gene}[2]<$pcut && $hyb2_par{$gene}<$pcut) # trans only
    {
      print OUT "\n$gene\t$hyb2{$gene}[0]\t$hyb2{$gene}[1]\t$hyb2{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb2_par{$gene}\ttrans";    
    }
    elsif($hyb2{$gene}[2]<$pcut && $par{$gene}[2]>=$pcut && $hyb2_par{$gene}<$pcut) # compensatory
    {
      print OUT "\n$gene\t$hyb2{$gene}[0]\t$hyb2{$gene}[1]\t$hyb2{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb2_par{$gene}\tcompensatory";    
    }
    elsif($hyb2{$gene}[2]>=$pcut && $par{$gene}[2]>=$pcut && $hyb2_par{$gene}>=$pcut) # conserved
    {
      print OUT "\n$gene\t$hyb2{$gene}[0]\t$hyb2{$gene}[1]\t$hyb2{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb2_par{$gene}\tconserved";    
    }
    else # ambiguous
    {
      print OUT "\n$gene\t$hyb2{$gene}[0]\t$hyb2{$gene}[1]\t$hyb2{$gene}[2]\t$par{$gene}[0]\t$par{$gene}[1]\t$par{$gene}[2]\t$hyb2_par{$gene}\tambiguous";    
    }
  }
  close OUT;
}

