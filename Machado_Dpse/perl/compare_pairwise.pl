#!/usr/bin/perl

########################################################################
# 
# 07/30/2011
#
# compare_pairwise.pl
#
# Purpose: compare FASTA features after they've been pairwise aligned
# 
# Input: two fasta files with some common features (chromosomes, exons, genes, etc.) 
#
# Output: lengths of common features, # that differ, etc.
#         
# Syntax: perl compare_pairwise.pl <FASTA1> <FASTA2>
#
########################################################################

use strict;
use warnings;

my$ref1 = $ARGV[0];
my$ref2 = $ARGV[1];
#my$out3 = $ARGV[2];

my$loc; my$seq; my$first; my$length1; my$length2; 
my$base1; my$base2; my$i; my$common; my$SNPs; my$perKB; my$j;
my$gaps1; my$gaps2; my$pos1; my$pos2;
my$out1; my$out2; my$out4; my$out5;
my$gene;

my@elements; my@seq1; my@seq2;

my%genome1; my%genome2; my%features;
my%diffs; my%coding; my%lengths; my%genes;

# read in reference sequences

$first = 1; $seq = "";
open(REF,"$ref1") or die "\nError opening $ref1\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome1{$loc} = $seq;} # the first time, don't save sequence in hash
    @elements = split(/\|/,$1);
    $loc = $elements[0];
    $features{$loc} = 1;
    $seq = "";
    $first = 0;
    #print "\nChromosome: $loc\n";
  }
  else
  {
    $seq = $seq.uc($_);
    #print "\nSequence: $seq\n";
  }
}
close REF;
$genome1{$loc} = $seq; # this needs to be done since last entry will not be stored in hash
#print "\n\nReference1 stored\n";

$first = 1; $seq = "";
open(REF,"$ref2") or die "\nError opening $ref2\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome2{$loc} = $seq;} # the first time, don't save sequence in hash
    @elements = split(/\|/,$1);
    $loc = $elements[0];
    $features{$loc} = 1;
    $seq = "";
    $first = 0;
    #print "\nChromosome: $loc\n";
  }
  else
  {
    $seq = $seq.uc($_);
    #print "\nSequence: $seq\n";
  }
}
close REF;
$genome2{$loc} = $seq; # this needs to be done since last entry will not be stored in hash
#print "\n\nReference2 stored\n";

if($ref1 =~ /^(\S+)\.(fasta|fa)/){$out1 = $1."\.SNPs.txt"; $out4 = $1."\_masked.fasta";}
if($ref2 =~ /^(\S+)\.(fasta|fa)/){$out2 = $1."\.SNPs.txt"; $out5 = $1."\_masked.fasta";}

#open(OUT1,"\> $out1") or die "Error opening $out1!";
#open(OUT2,"\> $out2") or die "Error opening $out2!";
#open(OUT3,"\> $out3") or die "Error opening $out3!";
#open(OUT4,"\> $out4") or die "Error opening $out4!";
#open(OUT5,"\> $out5") or die "Error opening $out5!";

# compare common features between genomes

#$SNPs = 0;

#print "locus\tmel\tsim\tfsa\n";

foreach $loc (keys %features)
{
  if($genome1{$loc} && $genome2{$loc})
  {
    $length1 = length($genome1{$loc});
    $length2 = length($genome2{$loc});

    #print "$loc\t$length1\t$length2\n";

    #print OUT4 "\>$loc\n";
    #print OUT5 "\>$loc\n"; 

    if($length1 == $length2)
    {
      if($loc =~ /(FBgn\d+)\_\S+/){$gene = $1;}
      #print "$gene\n"; exit;

      $genes{$gene} = 1;

      $lengths{$gene} += $length1; 

      #$common = 0;
      $SNPs = 0;
      $gaps1 = 0;
      $gaps2 = 0;

      for($i=0;$i<$length1;$i++)
      {
        $base1 = substr($genome1{$loc},$i,1);
        $base2 = substr($genome2{$loc},$i,1);
        $j = $i+1;

        if($base1 eq "-" || $base2 eq "-" || $base1 eq "N" || $base2 eq "N")
	{
          #print OUT4 "N";
          #print OUT5 "N";
        }
        else
	{
          #print OUT4 "$base1";
          #print OUT5 "$base2";
        }

        if($base1 eq "-"){$gaps1 += 1;}
        if($base2 eq "-"){$gaps2 += 1;}

        unless($base1 eq "-" || $base2 eq "-" || $base1 eq "N" || $base2 eq "N")
	{
          $coding{$gene} += 1;

          if($base1 ne $base2)
          {
            $diffs{$gene} += 1;
            $SNPs += 1;
            $pos1 = $j - $gaps1;
            $pos2 = $j - $gaps2;
            #print OUT1 "$loc\t$pos1\t$base1\t$base2\n";
            #print OUT2 "$loc\t$pos2\t$base2\t$base1\n";
            #print OUT3 "$loc\t$j\t$base1\t$base2\n";
            #print "$loc\t$pos1\t$pos2\t$j\n";
          }
	}
      }
      #$perKB = $SNPs/($common/1000);
      #print "$loc\t$length1\t$SNPs\n";
      #print OUT4 "\n";
      #print OUT5 "\n";
    }
    #else{print "$loc: \{$length1\}\t\{$length2\}\n";}
  }
}

print "gene\tlength\tcoding\tSNPs\n";

foreach $gene (keys %genes)
{
  if($diffs{$gene}){print "$gene\t$lengths{$gene}\t$coding{$gene}\t$diffs{$gene}\n";}
  else{print "$gene\t$lengths{$gene}\t$coding{$gene}\t0\n";}
}

#print "\nTotal coding sequence (without gaps): $totcod";
#print "\nTotal gapped sequence: $totgaps";
#print "\nNumber of SNPs: $SNPs\n\n";
