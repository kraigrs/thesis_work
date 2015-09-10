#!/usr/bin/perl

###################################################################################
#
# 07/05/2011
#
# imprinting_nullDist.pl
#
# Purpose: create a null distribution to compare original imprinting dataset to from 
#          both sliding window approaches
# 
# Input: which window approach to use, number of datasets to generate, how many genes to draw
#
# Output: data for each window, including number of genes (see output from respective window methods) 
#
# Fix: combined the two separate methods to view imprinting (gene clusters and genomic distance)
#
# Usage: perl imprinting_nullDist.pl <permutations> <# genes to draw> <const> <all> <inter> <folder> <window> <step> <k>
#        <permutations>    ==> number of permuted datasets to generate (~500?)
#        <# genes to draw> ==> number of genes to draw from original dataset
#        <const>           ==> list of constitutive genes
#        <all>             ==> list of all genes tested for imprinting 
#        <inter>           ==> all genes intersected with windows
#        <folder>          ==> where to put the data
#        <window>          ==> window size
#        <step>            ==> step to increment window
#        <k>               ==> size of total gene cluster
#
#  E.g. perl imprinting_nullDist.pl 2 130 ../McManus/constitutive_genes.txt ../mel_mel_data/zhr_z30_MB_cut20_genes.txt ../mel_mel_data/null_dist/all_genes_intersect.bed ../mel_mel_data/null_dist 500000 10000 50
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
  my$t0 = time;

  my$perms = $ARGV[0];
  my$draw = $ARGV[1];
  my$const = $ARGV[2];
  my$all = $ARGV[3];
  my$inter = $ARGV[4];
  my$prefix = $ARGV[5];

  my$wind = $ARGV[6];
  my$step = $ARGV[7];

  my$k = $ARGV[8];

  our($chr,$gene,$exon,$start,$stop,$genesRef,$impRef,$begin,$end,$key,$n,$oldGene,$max,$imps,$coord,$i,$window,$iter);
  our(@elements,@genes,@imps);
  our(%geneHash_gene,%geneHash_dist,%imprint,%all,%windows,%imp_wind,%dist_wind,%gene_wind);

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
      #print "gene=$gene\tstrand=$strand\n";
      if($chr)
      {
        if($chr ne "-")
	{
          #print "chr$chr\t$gene\t$start\t$stop\n";
          $geneHash_gene{$chr}{$start} = [$gene,$stop];
          $geneHash_dist{$chr}{$gene} = [$start,$stop];
          # stupid structure, but should allow easier sorting (by chromosome and position)
        }
        #else{print "$gene\n";}
      }
      #else{print "$gene\n";}
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

  foreach $chr (keys %geneHash_dist)
  {
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
      $coord = $start.",".$stop; 
      $windows{$chr}{$coord} = 0;
      $start += $step;
      $stop += $step;
    }
  }

  open(ALL,"$inter") or die "Can't open $inter for reading!";
  while(<ALL>) 
  {
    chomp;
    @elements = split(/\s+/,$_);
    $chr = $elements[4];
    $start = $elements[5];
    $stop = $elements[6];
    $coord = $start.",".$stop;

    $windows{$chr}{$coord} += 1;    
  } 
  close ALL;

  # start permutations

  $iter = 1;

  while($iter <= $perms)
  {
    system "R --slave --args $draw $all $prefix < ../Rscripts/imprinting_draw.r > ./trash";
  
    open(IMP,"$prefix/temp_imprint.txt") or die "Can't open $prefix/temp_imprint.txt for reading!";
    while(<IMP>) 
    {
      chomp;
      $imprint{$_} = 1;
    } 
    close IMP;
    
    #################
    #distance method#
    #################

    foreach $chr (keys %geneHash_dist)
    {
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
        $coord = $start.",".$stop; 
        $imp_wind{$chr}{$coord} = 0;
        $start += $step;
        $stop += $step;
      }
    }

    open(OUT,">$prefix\/imprinted_genes.bed") or die "Error writing to $prefix\/imprinted_genes.bed\n";
    foreach $gene (keys %imprint)
    {
      foreach $chr (keys %geneHash_dist)
      {
        if($geneHash_dist{$chr}{$gene}){print OUT "$chr\t$geneHash_dist{$chr}{$gene}[0]\t$geneHash_dist{$chr}{$gene}[1]\t$gene\n";}
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

      $imp_wind{$chr}{$coord} += 1;    
    } 
    close IMP;

    foreach $chr (keys %windows)
    {
      foreach $coord (keys %{$windows{$chr}})
      {
        $imps = $imp_wind{$chr}{$coord};
        $n = $windows{$chr}{$coord};      
        #print "$coord\n";
   
        $dist_wind{$chr}{$coord} = [$imps,$n];
      } 
    }

    # write output to file
    foreach $chr (keys %dist_wind)
    {    
      open(OUT1,">>$prefix\/chr$chr\/zhr_z30\.chr$chr\_wind$wind\_step$step\.imprinting.txt") or die "Error writing to $prefix\/chr$chr\/zhr_z30\.chr$chr\_wind$wind\_step$step\.imprinting.txt\n";
      if($iter == 1)
      {
        print OUT1 "coordinates";
        foreach $coord (sort {$a cmp $b} keys %{$dist_wind{$chr}})
        {
          print OUT1 "\t$coord";
        }

        print OUT1 "\ntotals";
        foreach $coord (sort {$a cmp $b} keys %{$dist_wind{$chr}})
        {
          print OUT1 "\t$dist_wind{$chr}{$coord}[1]";
        }

        print OUT1 "\nperm$iter";
        foreach $coord (sort {$a cmp $b} keys %{$dist_wind{$chr}})
        {
          print OUT1 "\t$dist_wind{$chr}{$coord}[0]";
        }
      }
      else
      {
        print OUT1 "\nperm$iter";
        foreach $coord (sort {$a cmp $b} keys %{$dist_wind{$chr}})
        {
          print OUT1 "\t$dist_wind{$chr}{$coord}[0]";
        }
      }
      close OUT1;
    }
    
    system "rm $prefix\/imprinted_genes.bed $prefix\/imprinted_genes_intersect.bed";
    undef %imp_wind; 
    #######end#######

    ################
    #cluster method#
    ################

    foreach $chr (keys %geneHash_gene)
    {
      $i = 0;
      foreach $start (sort {$a <=> $b} keys %{$geneHash_gene{$chr}})
      {
        #print "$i\n";
        if($i < $k)
        { 
          if($all{$geneHash_gene{$chr}{$start}[0]})
	  {
            push @genes,$geneHash_gene{$chr}{$start}[0];
            if($imprint{$geneHash_gene{$chr}{$start}[0]}){push @imps,$geneHash_gene{$chr}{$start}[0];}
            $i += 1;
            #print "$i\tchr=$chr\tgene=$geneHash_gene{$chr}{$start}[0]\n";
          }    
       
          if($i == $k)
          {
            $genesRef = \@genes;
            $impRef = \@imps;
            $n = &comp_array($genesRef,$impRef); #print "$n\n";
            $oldGene = shift @genes; #print "oldGene=$oldGene\tchr=$chr\n";
            if($imprint{$oldGene}){shift @imps;}
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
            $coord = $begin.",".$end;           

            $gene_wind{$chr}{$coord} = [$n,$i];
            $i -= 1;
          }
        }
      }
      undef @genes; undef @imps;
    }

    # write output to file
    foreach $chr (keys %gene_wind)
    {    
      open(OUT2,">>$prefix\/chr$chr\/zhr_z30\.chr$chr\_k$k.imprinting.txt") or die "Error writing to $prefix\/chr$chr\/zhr_z30\.chr$chr\_k$k.imprinting.txt\n";
      if($iter == 1)
      {
        print OUT2 "coordinates";
        foreach $coord (sort {$a cmp $b} keys %{$gene_wind{$chr}})
        {
          print OUT2 "\t$coord";
        }

        print OUT2 "\ntotals";
        foreach $coord (sort {$a cmp $b} keys %{$gene_wind{$chr}})
        {
          print OUT2 "\t$gene_wind{$chr}{$coord}[1]";
        }

        print OUT2 "\nperm$iter";
        foreach $coord (sort {$a cmp $b} keys %{$gene_wind{$chr}})
        {
          print OUT2 "\t$gene_wind{$chr}{$coord}[0]";
        }
      }
      else
      {
        print OUT2 "\nperm$iter";
        foreach $coord (sort {$a cmp $b} keys %{$gene_wind{$chr}})
        {
          print OUT2 "\t$gene_wind{$chr}{$coord}[0]";
        }
      }
      close OUT2;
    }
    #######end######

    $iter++;
    undef %imprint; 
  }
  printf ("\nTime elapsed: %d\n\n",time-$t0);
}

