#!/usr/bin/perl

######################################################################################
# 
# 10/27/2009
#
# geneSNPcoverage.pl
#
# Purpose: obtain summary information on the SNPs contained within a gene
# 
# Input: a text file with a list of genes and associated SNPs, with reference allele,
#        number of concordant (discordant) sequences, etc. This text file will be found 
#        based on whatever gene is typed by the user         
#
# Output: an output for the specified gene containing the SNP position, coverage,
#         allele counts, nucleotide counts
#
######################################################################################

use strict;
use warnings;

my$start = time;

my$data_dir = $ARGV[0];
my$gene = $ARGV[1];

my@list;
my$pos;
my$covseqs;
my$call;

my@alleles;
my$match;
my$As;
my$Cs;
my$Gs;
my$Ts;
my$mismatch;

open(OUT,"> $data_dir/$gene.summary.txt") or die "Error writing to $data_dir/$gene.summary.txt\n";
print OUT "pos\tcovseqs\tmatch\tmismatch\tA\tC\tG\tT";
close OUT;

open(GENE,"$data_dir/$gene.txt") or die "\nError opening $data_dir/$gene.txt\n";
while(<GENE>)
{
  @list = split(" ",$_);
  $pos = $list[1];
  $covseqs = $list[7];
  $call = $list[8];
  
  $match = 0;
  $As = 0;
  $Cs = 0;
  $Gs = 0;
  $Ts = 0; 
  $mismatch = 0;

  @alleles = split("",$call);
  foreach(@alleles)
  {
    if(/[\.,]/){$match+=1;}
    elsif(/[Aa]/){$As+=1;}
    elsif(/[Cc]/){$Cs+=1;}
    elsif(/[Gg]/){$Gs+=1;}
    elsif(/[Tt]/){$Ts+=1;}
  }
  
  $mismatch = $As + $Cs + $Gs + $Ts;
 
  open(OUT,">> $data_dir/$gene.summary.txt") or die "Error writing to $data_dir/$gene.summary.txt\n";
  print OUT "\n$pos\t$covseqs\t$match\t$mismatch\t$As\t$Cs\t$Gs\t$Ts";
  close OUT;
}
close GENE;
