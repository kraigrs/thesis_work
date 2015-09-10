#!/usr/bin/perl

###################################################################################
#
# 06/30/2011
#
# imprinting_dist.pl
#
# Purpose: take a list of imprinted genes and total genes tested for imprinting and look  
#          for enrichment based on chromosomal location
# 
# Input: intersected file, const exons, imprinted genes, all genes, prefix, perm 
#
# Output: number of genes in a given window 
#
# Usage: perl imprinting_dist.pl <wind> <step> <const> <imp> <intersect> <prefix> 
#        <wind>   ==> bp window 
#        <step>   ==> # bp to slide window
#        <const>  ==> list of constitutive genes!
#        <imp>    ==> list of imprinted genes
#        <intersect>   ==> intersected file of all genes 
#        <prefix> ==> what goes at beginning of file name
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
  my$wind = $ARGV[0];
  my$step = $ARGV[1];
  my$const = $ARGV[2];
  my$imp = $ARGV[3];
  my$inter = $ARGV[4];
  my$prefix = $ARGV[5];

  our($chr,$gene,$gene1,$gene2,$start,$stop,$max,$n,$imps,$coord);
  our(@elements);
  our(%geneHash,%imprint,%windows);

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
            $geneHash{$chr}{$gene} = [$start,$stop];
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
  print "\nRead $const\n";

  foreach $chr (keys %geneHash)
  {
    $max = 0;
    foreach $gene (keys %{$geneHash{$chr}})
    {
      if($max < $geneHash{$chr}{$gene}[1]){$max = $geneHash{$chr}{$gene}[1];} 
    }
    #print "max=$max\n";
    $start = 0; $stop = $wind;
    
    #create windows bedfile
    while($stop <= $max)
    {
      $coord = $start.",".$stop; 
      $windows{$chr}{$coord}{'total'} = 0;
      $windows{$chr}{$coord}{'imps'} = 0;
      $start += $step;
      $stop += $step;
    }
  }
  print "\nMade windows\n";
  
  open(ALL,"$inter") or die "Can't open $inter for reading!";
  while(<ALL>) 
  {
    chomp;
    @elements = split(/\s+/,$_);
    $chr = $elements[4];
    $start = $elements[5];
    $stop = $elements[6];
    $coord = $start.",".$stop;

    $windows{$chr}{$coord}{'total'} += 1;    
  } 
  close ALL;
  print "\nRead $inter\n";

  open(IMP,"$imp") or die "Can't open $imp for reading!";
  while(<IMP>) 
  {
    chomp;
    $imprint{$_} = 1;
  } 
  close IMP;
  print "\nRead $imp\n";

  #make imprinting list bed file to compare using intersectBed
  print "\nMaking $prefix\/imprinted_genes.bed\n";
  open(OUT,">$prefix\/imprinted_genes.bed") or die "Error writing to $prefix\/imprinted_genes.bed\n";
  foreach $gene (keys %imprint)
  {
    foreach $chr (keys %geneHash)
    {
      if($geneHash{$chr}{$gene}){print OUT "$chr\t$geneHash{$chr}{$gene}[0]\t$geneHash{$chr}{$gene}[1]\t$gene\n";}
    }
  }
  close OUT;

  system "intersectBed -a $prefix\/imprinted_genes.bed -b $prefix\/windows_dist.bed -wa -wb > $prefix\/imprinted_genes_intersect.bed";

  open(IMP,"$prefix\/imprinted_genes_intersect.bed") or die "Can't open $prefix\/imprinted_genes_intersect.bed for reading!";
  while(<IMP>) 
  {
    chomp;
    @elements = split(/\s+/,$_);
    $chr = $elements[4];
    $start = $elements[5];
    $stop = $elements[6];
    $coord = $start.",".$stop;

    $windows{$chr}{$coord}{'imps'} += 1;    
  } 
  close IMP;

  #foreach $chr (keys %windows)
  #{
  #  foreach $coord (keys %{$windows{$chr}})
  #  {
  #    print "chr$chr\tcoord=$coord\ttotal=$windows{$chr}{$coord}{'total'}\timps=$windows{$chr}{$coord}{'imps'}\n";
  #  }
  #}
  #exit;

  foreach $chr (keys %windows)
  {
    open(OUT,">$prefix\/chr$chr\/zhr_z30\.chr$chr\_wind$wind\_step$step\.imprinting.txt") or die "Error writing to $prefix\/chr$chr\/zhr_z30\.chr$chr\_wind$wind\_step$step\.imprinting.txt\n";
    print OUT "left_bound\tright_bound\ttotal\timprinted";
  
    foreach $coord (keys %{$windows{$chr}})
    {
      @elements = split(/,/,$coord);
      $imps = $windows{$chr}{$coord}{'imps'};
      $n = $windows{$chr}{$coord}{'total'};
      $start = $elements[0];
      $stop = $elements[1];      
      #print "$coord\n";
   
      print OUT "\n$start\t$stop\t$n\t$imps";
    } 
    close OUT;  
  }
  #system "rm $prefix\/imprinted_genes.bed $prefix\/imprinted_genes_intersect.bed";
}
