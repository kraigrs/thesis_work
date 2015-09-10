#!/usr/bin/perl

###################################################################################
#
# 06/29/2011
#
# imprinting_gene.pl
#
# Purpose: take a list of imprinted genes and total genes tested for imprinting and look  
#          for enrichment based on gene density
# 
# Input: list of imprinted genes, list of total number of genes tested, number of genes
#        to include in a window 
#
# Output: for a k-size gene window, the number of imprinted genes 
#
# Usage: perl imprinting_gene.pl <k> <const> <imp> <all> <prefix>
#        <k>       ==> number of genes 
#        <const>   ==> list of constitutive genes!
#        <imp>     ==> list of imprinted genes
#        <all>     ==> list of all genes tested for imprinting
#        <prefix>  ==> what goes at beginning of file name
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
  my$k = $ARGV[0];
  my$const = $ARGV[1];
  my$imp = $ARGV[2];
  my$all = $ARGV[3];
  my$prefix = $ARGV[4];

  our($chr,$gene,$exon,$start,$stop,$genesRef,$impRef,$begin,$end,$key,$n,$oldGene);
  my$i = 0;
  our(@elements,@genes,@imps);
  our(%geneHash,%imprint,%all);

  sub comp_array
  {
    my($genesRef,$impRef) = @_;
    my@a = @{$genesRef};
    my@b = @{$impRef};
    my$count = 0;
    my$item;
    my%counts;
    foreach $item (@a,@b){$counts{$item}++;}
    foreach $item (keys %counts)
    {
      if($counts{$item} == 2){$count += 1;}
      elsif($counts{$item} > 2){print "# matches exceeds 2! n = $counts{$item}\n"; exit;}
    }
    return $count;
  }

  # gene hash
  open(CONST,"$const") or die "Can't open $const for reading!";
  while(<CONST>) 
  {
    chomp;
    if($_ =~ /^#/){next;}
    else
    { 
      @elements = split(/\s+/,$_);
      $chr = $elements[1];
      $start = $elements[3];
      $stop = $elements[2];
      #print "$stop\n";
      $gene = $elements[0];
      if($chr)
      {
        if($chr ne "-")
	{
          #print "chr$chr\t$gene\t$start\t$stop\n";
          $geneHash{$chr}{$start} = [$gene,$stop];
          # stupid structure, but should allow easier sorting (by chromosome and position)
           
          if($start > $geneHash{$chr}{$start}[1]){print "Problem with gene=$gene: start=$start\tstop=$geneHash{$chr}{$start}[1]\n"; exit;}
        }
        #else{print "$gene\n";}
      }
      #else{print "$gene\n";}
    }
  }
  close CONST;
  #exit;

  open(IMP,"$imp") or die "Can't open $imp for reading!";
  while(<IMP>) 
  {
    chomp;
    $imprint{$_} = 1;
  } 
  close IMP;

  open(ALL,"$all") or die "Can't open $all for reading!";
  while(<ALL>) 
  {
    chomp;
    $all{$_} = 1;
  } 
  close ALL;

  foreach $chr (keys %geneHash)
  {
    open(OUT,">$prefix\/chr$chr\/zhr_z30\.chr$chr\_k$k\.imprinting.txt") or die "Error writing to $prefix\/chr$chr\/zhr_z30\.chr$chr\_k$k\.imprinting.txt\n";
    print OUT "left_bound\tright_bound\ttotal\timprint";
    foreach $start (sort {$a <=> $b} keys %{$geneHash{$chr}})
    {
      #print "$i\n";
      if($i < $k)
      { 
        if($all{$geneHash{$chr}{$start}[0]})
	{
          push @genes,$geneHash{$chr}{$start}[0];
          if($imprint{$geneHash{$chr}{$start}[0]}){push @imps,$geneHash{$chr}{$start}[0];}
          $i += 1;
          #print "$i\tchr=$chr\tgene=$geneHash{$chr}{$start}[0]\n";
        }  
       
        if($i == $k)
        {
          $genesRef = \@genes;
          $impRef = \@imps;
          $n = &comp_array($genesRef,$impRef); #print "$n\n";
          $oldGene = shift @genes; #print "oldGene=$oldGene\tchr=$chr\n";
          if($imprint{$oldGene}){shift @imps;}
          foreach $key (keys %{$geneHash{$chr}})
          {
            #print "$key\n";
            if($oldGene eq $geneHash{$chr}{$key}[0])
            {
              $begin = $key;
              #print "begin=$begin\n";
            }
          }
          $end = $geneHash{$chr}{$start}[1]; #print "endGene=$geneHash{$chr}{$start}[0]\tend=$end\n";
          if($begin && $end)
	  {
            print OUT "\n$begin\t$end\t$i\t$n";
          }
          else{die "Something is definitely going wrong here\n";}
          $i -= 1;
        }
      }
    }
    close OUT; $i = 0; undef @genes; undef @imps;
  }
}
