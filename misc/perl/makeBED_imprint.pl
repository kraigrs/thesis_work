#!/usr/bin/perl

###################################################################################
#
# 07/07/2011
#
# makeBED_imprint.pl
#
# Purpose: make BED files for the windows and all genes (to be run once before imprinting_dist.pl)
# 
# Input: window size, step size, const exons, all genes, prefix 
#
# Output: number of genes in a given window 
#
# Usage: perl makeBED_imprint.pl <wind> <step> <k> <const> <all> <prefix>
#        <wind>   ==> bp window 
#        <step>   ==> # bp to slide window
#        <k>      ==> number of genes per cluster
#        <const>  ==> list of constitutive genes!
#        <all>    ==> list of imprinted genes
#        <prefix> ==> what goes at beginning of file name
#
# E.g. perl makeBED_imprint.pl 500000 10000 50 ../McManus/constitutive_genes.txt ../mel_mel_data/imprinted_genes_112211.txt ../mel_mel_data/zhr_z30_genes.txt ../mel_mel_data/imprinting
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
  my$wind = $ARGV[0];
  my$step = $ARGV[1];
  my$k = $ARGV[2];
  my$const = $ARGV[3];
  my$imps = $ARGV[4];
  my$all = $ARGV[5];  
  my$prefix = $ARGV[6];

  our($chr,$gene,$gene1,$gene2,$start,$stop,$max,$oldGene,$key,$begin,$end);
  my$i = 0;
  our(@elements,@genes);
  our(%geneHash_dist,%geneHash_gene,%all);

  # gene hash
  open(CONST,"$const") or die "Can't open $const for reading!";
  while(<CONST>) 
  {
    chomp;
    if($_ =~ /^#/){next;}
    else
    { 
      @elements = split(/\s+/,$_);
      $chr = $elements[2];
      $start = $elements[4]; # subtract 1 to make it 0-based for BED format
      $stop = $elements[3];
      #print "$stop\n";
      $gene1 = $elements[0];
      $gene2 = $elements[1];
      #print "gene=$gene\tstrand=$strand\n";
      unless(!$gene2)
      {
        if($gene2 ne "-")
        {
          $gene = $gene1;
          if($chr ne "-")
	  {
            #print "chr$chr\t$gene\t$start\t$stop\n";
            $start = $start - 1;  # subtract 1 to make it 0-based for BED format
            $geneHash_gene{$chr}{$start} = [$gene,$stop];
            $geneHash_dist{$chr}{$gene} = [$start,$stop];
            # stupid structure, but should allow easier sorting (by chromosome and position)
          }
          #else{print "$gene\n";}
        }
        #else{print "$gene\n";}
      }
    }
  }
  close CONST;
  #exit;

  open(ALL,"$all") or die "Can't open $all for reading!";
  while(<ALL>) 
  {
    chomp;
    $all{$_} = 1;
  } 
  close ALL;

  #make imprinting and all list bed files to compare using intersectBed

  open(OUT,">$prefix\/all_genes.bed") or die "Error writing to $prefix\/all_genes.bed\n";
  foreach $gene (keys %all)
  {
    foreach $chr (keys %geneHash_dist)
    {
      if($geneHash_dist{$chr}{$gene}){print OUT "$chr\t$geneHash_dist{$chr}{$gene}[0]\t$geneHash_dist{$chr}{$gene}[1]\t$gene\n";}
    }
  }
  close OUT;
  #exit;

  open(OUT,">$prefix\/windows_dist.bed") or die "Error writing to $prefix\/windows_dist.bed\n";  
  foreach $chr (keys %geneHash_dist)
  {
    open(TEMP,">$prefix\/chr$chr\/chr$chr\_windows_dist.bed") or die "Error writing to $prefix\/chr$chr\/chr$chr\_windows_dist.bed\n";
    $max = 0;
    foreach $gene (keys %{$geneHash_dist{$chr}})
    {
      if($max < $geneHash_dist{$chr}{$gene}[1]){$max = $geneHash_dist{$chr}{$gene}[1];} 
    }
    #print "max=$max\n";
    $start = 0; $stop = $wind;
    
    #create windows bedfile
     
    while($stop <= $max)
    { 
      print TEMP "$chr\t$start\t$stop\n";
      print OUT "$chr\t$start\t$stop\n";
      $start += $step;
      $stop += $step;
    }
    close TEMP;
  }
  close OUT;

  system "intersectBed -a $prefix\/all_genes.bed -b $prefix\/windows_dist.bed -wa -wb > $prefix\/all_genes_intersect.bed";

  open(TEMP,">$prefix\/windows_gene.bed") or die "Error writing to $prefix\/windows_gene.bed\n"; 
  foreach $chr (keys %geneHash_gene)
  {
    #print "Just past temp\n";
    open(OUT,">$prefix\/chr$chr\/chr$chr\_windows_gene.bed") or die "Error writing to $prefix\/chr$chr\/chr$chr\_windows_gene.bed\n";
    foreach $start (sort {$a <=> $b} keys %{$geneHash_gene{$chr}})
    {
      #print "Just past out\n";
      if($i < $k)
      { 
        if($all{$geneHash_gene{$chr}{$start}[0]})
	{
          #print "Incrementing gene list\n";
          push @genes,$geneHash_gene{$chr}{$start}[0];
          $i += 1;
          #print "i=$i\tk=$k\n";
        }  
      }
      #print "i=$i\tk=$k\n";
      if($i == $k)
      {
        $oldGene = shift @genes; #print "oldGene=$oldGene\tchr=$chr\n";
        #print "Gene cluster is saturated\n";

        foreach $key (keys %{$geneHash_gene{$chr}})
        {
          #print "$key\n";
          if($oldGene eq $geneHash_gene{$chr}{$key}[0])
          {
            $begin = $key;
            #print "begin=$begin\n";
          }
        }
        $end = $geneHash_gene{$chr}{$start}[1]; #print "endGene=$geneHash_gene{$chr}{$start}[0]\tend=$end\n";

        if($begin && $end)
        {
          #print "Should be printing some stuff";
          print OUT "$chr\t$begin\t$end\n";
          print TEMP "$chr\t$begin\t$end\n";
        }
        else{die "Something is definitely going wrong here\n";}
        $i -= 1;
      }
    }
    close OUT; $i = 0; undef @genes;
  }
  close TEMP;
}
