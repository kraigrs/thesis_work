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
# Syntax: perl compare_pairwise.pl <FASTA1> <FASTA2> <output>
#
########################################################################

use strict;
use warnings;

my$ref1 = $ARGV[0];
my$ref2 = $ARGV[1];
#my$out3 = $ARGV[2];

my$loc; my$seq; my$first; my$length1; my$length2; my$base1; my$base2; my$i; my$common; my$SNPs; my$perKB; my$j;
my$gaps1; my$gaps2; my$pos1; my$pos2;
my$out1; my$out2; my$out4; my$out5;
my@elements; my@seq1; my@seq2;
my%genome1; my%genome2; my%features;

# read in reference sequences

$first = 1; $seq = "";
open(REF,"$ref1") or die "\nError opening $ref1\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome1{$loc} = $seq;} # the first time, don't save sequence in hash
    #@elements = split(/\|/,$1);
    $loc = $1;
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
#exit;
$first = 1; $seq = "";
open(REF,"$ref2") or die "\nError opening $ref2\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)$/)
  {
    unless($first == 1){$genome2{$loc} = $seq;} # the first time, don't save sequence in hash
    #@elements = split(/\|/,$1);
    $loc = $1;
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

if($ref1 =~ /^(\S+)\.(fasta|fa)/){$out1 = $1."\.SNPs.txt"; $out4 = $1."\_masked.fa";}
if($ref2 =~ /^(\S+)\.(fasta|fa)/){$out2 = $1."\.SNPs.txt"; $out5 = $1."\_masked.fa";}

#open(OUT1,"\> $out1") or die "Error opening $out1!";
#open(OUT2,"\> $out2") or die "Error opening $out2!";
#open(OUT3,"\> $out3") or die "Error opening $out3!";
#open(OUT4,"\> $out4") or die "Error opening $out4!";
#open(OUT5,"\> $out5") or die "Error opening $out5!";

# compare common features between genomes

print "locus\tSNPs\tcommon\n";

#$SNPs = 0;
foreach $loc (keys %features)
{
  if($genome1{$loc} && $genome2{$loc})
  {
    #print "here\n";
    $length1 = length($genome1{$loc});
    $length2 = length($genome2{$loc});
    #if($length1 != $length2){print "$loc\n";}

    $common = 0;
    $SNPs = 0;

    $gaps1 = 0;
    $gaps2 = 0;

    #print "$loc\t$length1\t$length2\n";

    #print OUT4 "\>$loc\n";
    #print OUT5 "\>$loc\n"; 

    if($length1 == $length2)
    {
      for($i=0;$i<$length1;$i++)
      {
        $base1 = uc(substr($genome1{$loc},$i,1));
        $base2 = uc(substr($genome2{$loc},$i,1));
        $j = $i+1;

        if($base1 =~ /[^ACGTacgt]/ || $base2 =~ /[^ACGTacgt]/)
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

        unless($base1 =~ /[^ACGTacgt]/ || $base2 =~ /[^ACGTacgt]/)
	{
          #print "Base1: $base1\tBase2: $base2\n";
          $common += 1;
          if($base1 ne $base2)
          {
            $SNPs += 1;
            $pos1 = $j - $gaps1;
            $pos2 = $j - $gaps2;
            #print OUT1 "$loc\t$pos1\t$base1\t$base2\n";
            #print OUT2 "$loc\t$pos2\t$base2\t$base1\n";
            #print OUT3 "$loc\t$pos1\t$pos2\t$j\t$base1\t$base2\n";
            #print OUT1 "$loc\_$SNPs\t$pos1\t$base1\t$base2\n";
            #print OUT2 "$loc\_$SNPs\t$pos2\t$base2\t$base1\n";
            #print OUT3 "$loc\_$SNPs\t$j\t$base1\t$base2\n";
          }
	}
      }
      #$perKB = $SNPs/($common/1000);
      #print "$loc\t$length1\t$SNPs\n";
      #print OUT4 "\n";
      #print OUT5 "\n";

      print "$loc\t$SNPs\t$common\n";
    }
    else{print "$loc\n";}
  }
}

#print "$SNPs\n";
