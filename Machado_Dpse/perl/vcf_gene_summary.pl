#!/usr/bin/perl

########################################################################
# 
# 03/21/2014
#
# vcf_gene_summary.pl
#
# Purpose: for each gene, summarize variants that passed and failed filters, and heterozygosity
# 
# Input: VCF file 
#
# Output: list of genes and properties of the variants called
#         
# Syntax: perl vcf_gene_summary.pl <VCF file> 
#
########################################################################

use strict;
use warnings;
use List::Util qw(min max);

my$bed = $ARGV[0];
my$vcf = $ARGV[1];
#my$output = $ARGV[2];

my$gene; my$line; my$locus; 
my@elements;
my%gene_list; 

open(BED,"$bed") or die "Can't open $bed for reading!";
while(<BED>) 
{
  chomp;
  unless(/^#/)
  { 
    @elements = split(/\s/,$_);
    $gene = $elements[3];
    $gene_list{$gene}{PASS}{HET} = 0;
    $gene_list{$gene}{PASS}{HOM} = 0;
    $gene_list{$gene}{FAIL}{HET} = 0;
    $gene_list{$gene}{FAIL}{HOM} = 0;
  }
}
close BED;

open(VCF,"$vcf") or die "\nError opening $vcf\n";
while(<VCF>)
{
  chomp;
  $line = $_;
  unless(/^#/)
  {
    @elements = split(/\s+/,$_);
    $locus = $elements[0];
    if($locus =~ /(FBgn\d+)\_/){$gene = $1;}

    if($elements[6] =~ /PASS/)
    {
      if($elements[7] =~ /\S*(AF=(1\.00|0\.500|0\.500,0\.500));/)
      {
        if($2 eq "1\.00")
        {
          $gene_list{$gene}{PASS}{HOM} += 1;
        }
        else
        {
          $gene_list{$gene}{PASS}{HET} += 1;
        }
      } 
    }
    elsif($elements[6] =~ /FAIL/)
    {      
      if($elements[7] =~ /\S*(AF=(1\.00|0\.500|0\.500,0\.500));/)
      {
        if($2 eq "1\.00")
        {
          $gene_list{$gene}{FAIL}{HOM} += 1;
        }
        else
        {
          $gene_list{$gene}{FAIL}{HET} += 1;
        }
      }
    }
  }
}
close VCF;

print "gene\tpass_hom\tpass_het\tfail_hom\tfail_het\n";
foreach $gene (keys %gene_list)
{
  print "$gene\t$gene_list{$gene}{PASS}{HOM}\t$gene_list{$gene}{PASS}{HET}";
  print "\t$gene_list{$gene}{FAIL}{HOM}\t$gene_list{$gene}{FAIL}{HET}\n";
}
