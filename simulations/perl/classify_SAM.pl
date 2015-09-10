#!/usr/bin/perl

###################################################################################
#
# 04/06/2012
#
# classify_SAM.pl
#
# Purpose: using existing SNP information, assign reads to an allele
# 
# Input: SAM file that has been aligned to a single reference allowing mismatches 
#
# Output: a list of genes and their respective match and mismatch counts, along with a 
#         summary of all the reads 
#
# E.g. perl classify_SAM.pl 
#
###################################################################################

use strict;
use warnings;
use Switch;

my$ref = $ARGV[0];
my$const = $ARGV[1];
my$species1 = $ARGV[2];
my$species2 = $ARGV[3];
my$SNPs = $ARGV[4];
my$sam = $ARGV[5];
my$prefix = $ARGV[6];

our($chr,$pos,$base,$name,$refct,$altct,$start,$stop,$flag,$n,$amb,$var,$gene,$exon,$gene_exon,$test1,$test2,$test3);
our(@elements,@read,@meta);
our(%SNPhash,%refHash,%geneHash,%exonHash);

my$found = 0; my$seq = "";

#print "\nBuilding reference hash...\n";

open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  #print "Here!\n";
  if(/^>(chr\S+)$/ && $found == 1)
  {
    $refHash{$chr} = $seq;

    $seq = "";
    $chr = $1;
  }
  elsif(/^>(chr\S+)$/)
  {
    $seq = "";
    $chr = $1;
    $found = 1;
  }
  elsif($found == 1)
  {
    $seq = $seq.uc($_);
  }
}
$refHash{$chr} = $seq;
close REF;

#foreach(keys %refHash)
#{
#  print "\nChromosome: $_\n";
#  $seq = substr($refHash{$_},0,50);
#  print "First 50 bases: $seq\n";
#}
#print "\n";
#exit;

# gene and exon hashes 
open(CONST,"$const") or die "Can't open $const for reading!";
LINE:while(<CONST>) 
{
  chomp;
  if($_ =~ /^track.+/){next LINE;}
  else
  { 
    @elements = split(/\s/,$_);
    $gene = $elements[3];
    #$exon = $elements[1] . "," . $elements[2];
    $exon = $elements[1] . "_" . $elements[2];

    $geneHash{$gene}{$species1} = 0; 
    $exonHash{$gene}{$exon}{$species1} = 0;

    $geneHash{$gene}{$species2} = 0; 
    $exonHash{$gene}{$exon}{$species2} = 0;

    $geneHash{$gene}{b} = 0; 
    $exonHash{$gene}{$exon}{b} = 0;
  }
}
close CONST;

#print "\nBuilding SNP hash...\n";
open(SNPS,"$SNPs") or die "\nError opening $SNPs\n";
while(<SNPS>)
{
  chomp;
  @elements = split("\t",$_);
  $chr = $elements[0];
  $start = $elements[1];
  $stop = $elements[2];
  $base = $elements[3];
  $gene = $elements[7];
  $exon = $elements[5]."_".$elements[6];
  $gene_exon = $gene."_".$exon;
  $pos = $start + 1;
  $SNPhash{$gene_exon}{$pos} = [$chr,$base];
}
close SNPS;

open(READS,">$prefix\.reads.txt") or die "Error writing to $prefix\.reads.txt\n";
print READS "read\tgene_exon\tcall\n";

#print "\nReading SAM file...\n";
open(SAM,"$sam") or die "\nError opening $sam\n";
while (<SAM>) 
{
  chomp;
  $refct = 0; $altct = 0;
 
  @elements = split(/\s/,$_);

  $gene = $elements[0];
  $exon = $elements[1]."_".$elements[2];
  $gene_exon = $gene."_".$exon;
  @meta = split(/\@/,$elements[3]);

  $name = $meta[0];
  $seq = $meta[1];
  $start = $meta[2] + 1;
  $stop = $meta[3];
  @read = split("",$seq);

  foreach $pos (keys %{$SNPhash{$gene_exon}}) 
  {
    if($start <= $pos && $pos <= $stop)
    {
      $chr = $SNPhash{$gene_exon}{$pos}[0];

      if($read[$pos-$start] eq substr($refHash{$chr},$pos-1,1)){$refct += 1;}
      elsif($read[$pos-$start] eq $SNPhash{$gene_exon}{$pos}[1]){$altct += 1;}
      #elsif($read[$pos-$start] eq &altAllele(substr($refHash{$chr},$pos-1,1),$SNPhash{$gene_exon}{$pos}[1]))
      #{
        #$altct += 1;
        #$test1 = $SNPhash{$gene_exon}{$pos}[1];
        #$test2 = substr($refHash{$chr},$pos-1,1);
        #$test3 = &altAllele(substr($refHash{$chr},$pos-1,1),$SNPhash{$gene_exon}{$pos}[1]);

        #print "Read base: $read[$pos-$start]\n";
        #print "SNP: $test1\n";
        #print "Reference allele: $test2\n";
        #print "Alternative allele: $test3\n";
        #exit;
      #}
    }
  }
    
  if($refct > 0 && $altct == 0)
  {
    $flag = $species1;
    $geneHash{$gene}{$species1} += 1; 
    $exonHash{$gene}{$exon}{$species1} += 1;
  }
  elsif($refct == 0 && $altct > 0)
  {
    $flag = $species2;
    $geneHash{$gene}{$species2} += 1; 
    $exonHash{$gene}{$exon}{$species2} += 1;
  }
  elsif($refct == 0 && $altct == 0)
  {
    $flag = "Both";
    $geneHash{$gene}{b} += 1; 
    $exonHash{$gene}{$exon}{b} += 1;
  }
  else{$flag = "error";}

  print READS "$name\t$gene_exon\t$flag\n";
  #exit;
}
close READS;

open(EXONS,">$prefix\.exons.txt") or die "Error writing to $prefix\.exons.txt\n";
print EXONS "gene_exon\t$species1\t$species2\tBoth\tAdj$species1\tAdj$species2";

my$madj = 0;
my$sadj = 0;

foreach my$k1 (keys %exonHash)
{
  foreach my$k2 (keys %{$exonHash{$k1}})
  {
    if($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2} == 0)
    {
      print EXONS "\n$k1\_$k2\t$exonHash{$k1}{$k2}{$species1}\t$exonHash{$k1}{$k2}{$species2}\t$exonHash{$k1}{$k2}{b}\t0\t0";
    }
    else
    {
      $madj = $exonHash{$k1}{$k2}{$species1} + ($exonHash{$k1}{$k2}{$species1} / ($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2})) * $exonHash{$k1}{$k2}{b};  
      $sadj = $exonHash{$k1}{$k2}{$species2} + ($exonHash{$k1}{$k2}{$species2} / ($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2})) * $exonHash{$k1}{$k2}{b};
      print EXONS "\n$k1\_$k2\t$exonHash{$k1}{$k2}{$species1}\t$exonHash{$k1}{$k2}{$species2}\t$exonHash{$k1}{$k2}{b}\t$madj\t$sadj";
    }
  }
}
close EXONS;

open(GENES,">$prefix\.genes.txt") or die "Error writing to $prefix\.genes.txt\n";
print GENES "gene\t$species1\t$species2\tBoth\tAdj_$species1\tAdj_$species2";

$madj = 0;
$sadj = 0;

foreach(keys %geneHash)
{
  if($geneHash{$_}{$species1} + $geneHash{$_}{$species2} == 0)
  {
    print GENES "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}\t$geneHash{$_}{b}\t0\t0";
  }
  else
  {
    $madj = $geneHash{$_}{$species1} + ($geneHash{$_}{$species1} / ($geneHash{$_}{$species1} + $geneHash{$_}{$species2})) * $geneHash{$_}{b};
    $sadj = $geneHash{$_}{$species2} + ($geneHash{$_}{$species2} / ($geneHash{$_}{$species1} + $geneHash{$_}{$species2})) * $geneHash{$_}{b};
    print GENES "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}\t$geneHash{$_}{b}\t$madj\t$sadj";
  }
}
close GENES;
  
sub altAllele
{
  ($n,$amb) = @_;
  switch ($amb)
  {
    case "R"
    {
      if($n eq "A"){$var = "G";}
      elsif($n eq "G"){$var = "A";}
      else{die "Something is seriously wrong here!\n";}
    }
    case "Y"
    {
      if($n eq "C"){$var = "T";}
      elsif($n eq "T"){$var = "C";}
      else{die "Something is seriously wrong here!\n";}
    }
    case "M"
    {
      if($n eq "A"){$var = "C";}
      elsif($n eq "C"){$var = "A";}
      else{die "Something is seriously wrong here!\n";}
    }
    case "K" 
    {
      if($n eq "G"){$var = "T";}
      elsif($n eq "T"){$var = "G";}
      else{die "Something is seriously wrong here!\n";}
    }
    case "S"
    {
      if($n eq "G"){$var = "C";}
      elsif($n eq "C"){$var = "G";}
      else{die "Something is seriously wrong here!\n";}
    }
    case "W"
    {
      if($n eq "A"){$var = "T";}
      elsif($n eq "T"){$var = "A";}
      else{die "Something is seriously wrong here!\n";}
    }
    else{die "Something is seriously wrong here!\n";}
  }
  return $var;
}
