#!/usr/bin/perl

######################################################################################
# 
# 12/09/2013
#
# convertVCF.pl
#
# Purpose: read in pileup from SAMtools to call SNPs (or mutations)
# 
# Input: SAMtools alignment pileup and VCF file         
#
# Output: based on certain criteria, return highest confident variants
#
######################################################################################

use strict;
use warnings;

my$fasta = $ARGV[0];
my$regions = $ARGV[1];
my$VCF = $ARGV[2];
my$output = $ARGV[3];

my$seq; my$first; my$chr; my$pos; my$ref; my$alt; my$gene; my$line; my$start;
my$qual; my$filter; my$info; my$new_pos;
my@elements;
my%genome; my%vcf_gene;

$seq=""; $first=1;
open(REF,"$fasta") or die "\nError opening $fasta\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)\s*/)
  {
    unless($first == 1){$genome{$chr} = $seq;} # the first time, don't save sequence in hash
    $chr = $1;
    $seq = "";
    $first = 0;
  }
  else
  {
    $seq = $seq.uc($_);
  }
}
$genome{$chr} = $seq; # this needs to be done since last entry will not be stored in hash
close REF;

# create intersection of VCF and regions file

system "perl vcf2bed.pl $VCF > $VCF.bed";

system "intersectBed -a $VCF.bed -b $regions -wa -wb -f 1 | ".
         "awk '{print \$1 \"\t\" \$2 \"\t\" \$3 \"\t\" \$5 \"\t\" \$6 \"\t\" \$7}' ".
         "> $VCF\_regions.bed";

open(REG,"$VCF\_regions.bed") or die "\nError opening $VCF\_regions.bed\n";
while(<REG>)    
{
  chomp;
  @elements = split("\t",$_);
  $chr = $elements[0];
  $pos = $elements[2];
  $start = $elements[3];
  $gene = $elements[5];
  $vcf_gene{$chr}{$pos} = [$gene,$start];
}
close REG;

open(OUT,"> $output") or die "Error writing to $output\n";

open(VCF,"$VCF") or die "\nError opening $VCF\n";
while(<VCF>)
{
  chomp;
  $line = $_;
  if(/^#/){print OUT "$line\n";}
  else
  {
    @elements = split("\t",$line);
    $chr = $elements[0];
    $pos = $elements[1];
    $ref = $elements[3];
    $alt = $elements[4];
    $qual = $elements[5];
    $filter = $elements[6];
    $info = $elements[7];

    if($vcf_gene{$chr}{$pos})
    {
      $gene = $vcf_gene{$chr}{$pos}[0];
      $start = $vcf_gene{$chr}{$pos}[1];
      $new_pos = $pos-$start;
      print OUT "$gene\t$new_pos\t\.\t$ref\t$alt\t$qual\t$filter\t$info\n";
    }    
  }
}
close VCF;

close OUT;
