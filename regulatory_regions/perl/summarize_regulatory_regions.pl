#!/usr/bin/perl

###################################################################################
#
# 03/06/2013
#
# summarize_regulatory_regions.pl
#
###################################################################################

use strict;
use warnings;

my$regions = $ARGV[0];
my$SNPs = $ARGV[1];

my$gene; my$loc; my$gene1; my$gene2; my$line; my$key1; my$key2; my$tot_length_intron; my$tot_num_intron;
my$tot_length_inter; my$tot_num_inter;
my@elements; my@meta; my@first; my@second;
my%all; my%regs; my%intron; my%inter;

open(REGS,"$regions") or die "Can't open $regions for reading!";
while(<REGS>) 
{
  chomp;
  $line = $_;
  @elements = split(/\s+/,$line);  # splits up a string based on a delimiter, in this case, white space
  $loc = $elements[3];
  $regs{$loc} = $elements[2] - $elements[1]; # save the length of each region
}
close REGS;

open(SNPS,"$SNPs") or die "Can't open $SNPs for reading!";
while(<SNPS>) 
{
  chomp;
  $line = $_;
  @elements = split(/\s+/,$line);  # splits up a string based on a delimiter, in this case, white space
  $loc = $elements[0];

  if($loc =~ /intergenic/)
  {
    @meta = split(/\_/,$loc);
    #foreach(@meta){print "$_\n";} exit;

    if(scalar(@meta) == 2)              # scalar counts the number of elements in an array
    {
      #print "You shouldnt be here!!!"; exit;
      if($meta[0] =~ /intergenic/ && $meta[1] =~ /CG\d+/)
      {
        #$inter{$meta[1]}{$loc} += 1;
        $inter{$meta[1]} += 1;

        $all{$meta[1]} = 1;
      }
      elsif($meta[0] =~ /CG\d+/ && $meta[1] =~ /intergenic/)
      {
        #$inter{$meta[0]}{$loc} += 1;
        $inter{$meta[0]} += 1;

        $all{$meta[0]} = 1;
      }
    }
    elsif(scalar(@meta) == 3)
    {
      #print "You are here!!!"; exit;
      if($meta[0] =~ /CG\d+/ && $meta[2] =~ /CG\d+/)
      {
        $all{$meta[0]} = 1; $all{$meta[2]} = 1;

        #print "$meta[0]\t$meta[2]\n"; exit;
        #$inter{$meta[0]}{$loc} += 1;
        #$inter{$meta[2]}{$loc} += 1;
        $inter{$meta[0]} += 1;
        $inter{$meta[2]} += 1;
        #print "Locus: $loc\tGene: $meta[0]\tSNPs: $inter{$meta[0]}{$loc}\n";
        #print "Locus: $loc\tGene: $meta[2]\tSNPs: $inter{$meta[2]}{$loc}\n"; exit;
      }
      elsif($meta[0] =~ /CG\d+/)
      {
        $all{$meta[0]} = 1;
        #$inter{$meta[0]}{$loc} += 1;
        $inter{$meta[0]} += 1;
      }
      elsif($meta[2] =~ /CG\d+/)
      {
        $all{$meta[2]} = 1;
        #$inter{$meta[2]}{$loc} += 1;
        $inter{$meta[2]} += 1;
      }
    }
  }
  elsif($loc =~ /intron/)
  {
    @meta = split(/\_/,$loc);
    @first = split(/\:/,$meta[1]);
    @second = split(/\:/,$meta[2]);
    $gene1 = $first[0];
    $gene2 = $second[0];
    
    if($gene1 eq $gene2)
    {
      $all{$gene1} = 1;
      #$intron{$gene1}{$loc} += 1;
      $intron{$gene1} += 1;
    }
  }  
}
close SNPS;

#print "region\tgene\tsites\tlength\n";
print "gene\tintronic_SNPs\tintergenic_SNPs\n";


foreach $gene (keys %all)
{
  #print "\{$gene\}\n";
  #$tot_length_inter = 0; $tot_num_inter = 0;
  #$tot_length_intron = 0; $tot_num_intron = 0;

  if($inter{$gene} && $intron{$gene})
  {
    print "$gene\t$intron{$gene}\t$inter{$gene}\n";
  }
  elsif($inter{$gene})
  {
    print "$gene\t0\t$inter{$gene}\n";

    #foreach $key1 (keys %{ $inter{$gene} })
    #{
    #  if($regs{$key1})
    #  {
    #    #print "$key1\n";
    #    $tot_length_inter += $regs{$key1}; $tot_num_inter += 1;
    #    #print "$key1\t$gene\t$inter{$gene}{$key1}\t$regs{$key1}\n";
    #  }
    #}
  }
  elsif($intron{$gene})
  {
    print "$gene\t$intron{$gene}\t0\n";

    #foreach $key2 (keys %{ $intron{$gene} })
    #{
    #  if($regs{$key2})
    #  {
    #    #print "$key2\n";
    #    $tot_length_intron += $regs{$key2}; $tot_num_intron += 1;
    #    #print "$key2\t$gene\t$intron{$gene}{$key2}\t$regs{$key2}\n";
    #  }
    #}
  }
  else{print "$gene\t0\t0\n";}

  #print "$gene\t$tot_num_intron\t$tot_num_inter\n";
}
