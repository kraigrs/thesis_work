#!/usr/bin/perl

###################################################################################
#
# 02/12/2009
#
# classify2.pl
#
# Purpose: (adapted from classify1.pl) classify reads (single and paired-end) as matching
#          or mismatching the combined reference genome. The difference here is not using 
#          SNP information, but rather simply aligning to a particular species in the reference
#          i.e. the reference contains D. mel. and D. sec. references, whichever a read aligns to,
#          it belongs to that species
# 
# Input: two separate alignment files (mates 1 and mates 2) 
#
# Output: a list of genes and their respective match and mismatch counts, along with a 
#         summary of all the reads 
#
# E.g. perl classify2.pl ../names.txt Dmel Dsec ../mel_sec_data/Hyb_mate1.dmelXdsec-gene.sam ../mel_sec_data/Hyb_mate2.dmelXdsec-gene.sam ../mel_sec_data/Hyb.dmelXdsec-gene.reads.txt ../mel_sec_data/Hyb.dmelXdsec-gene.genes.txt &
#
###################################################################################

use strict;
#use warnings;

my@elements;
my%sech2mel;

my$annotation = $ARGV[0];

open(NAMES,"$annotation") or die "Can't open $annotation for reading!";
while(<NAMES>)
{
  if(/^ABREV.+/){next;}
  else
  {
    @elements = split(/\s+/,$_); 
    unless(scalar@elements < 4)
    {
      #print "$elements[3]\t$elements[1]\n";
      $sech2mel{$elements[3]} = $elements[1];
    }
  }
}
close NAMES;

my$line1; my$line2;
my$read1; my$read2;
my$align1; my$align2;
my$gene1; my$gene2;
my$spec1; my$spec2;
my$seq1; my$seq2;
my@elements1; my@elements2;
my%geneHash;

my$species1 = $ARGV[1]; # Dmel
my$species2 = $ARGV[2]; # Dsec
my$SAMfile1 = $ARGV[3];
my$SAMfile2 = $ARGV[4];
my$read_out = $ARGV[5]; 

open(OUT,">> $read_out") or die "Error writing to $read_out\n";
print OUT "\nread\tseq1\tseq2\tgene1\tspecies1\tgene2\tspecies2\tcall";
close OUT;

open(SAM1,"$SAMfile1") or die "Can't open $SAMfile1 for reading!";
open(SAM2,"$SAMfile2") or die "Can't open $SAMfile2 for reading!";
while (<SAM1>) 
{
  chomp;
  $line1 = $_;
  $line2 = <SAM2>; 
  
  if($line1 =~ /^\@.+/){next;}
  elsif($line1 =~ /^HWI.+\/1\s+/)
  { 
    @elements1 = split(/\s/,$line1);
    if($elements1[0] =~ /^(HWI.+)\/\d{1}/){$read1 = $1;}
    $align1 = $elements1[2];
    if($align1 =~ /^(.+)\/(.+)$/)
    {
      $gene1 = $1; $spec1 = $2;
      #print "\nmate1\t$gene1\t$spec1";
    }
    else{$gene1 = "*"; $spec1 = "*";}
    $seq1 = $elements1[9];
    #print "$seq1\t";
  
    @elements2 = split(/\s/,$line2);
    if($elements2[0] =~ /^(HWI.+)\/\d{1}/){$read2 = $1;}
    $align2 = $elements2[2];
    if($align2 =~ /^(.+)\/(.+)$/)
    {
      $gene2 = $1; $spec2 = $2;
      #print "\nmate2\t$gene2\t$spec2";
    }
    else{$gene2 = "*"; $spec2 = "*";}
    $seq2 = $elements2[9];
    #print "$seq2\n";

    if(($read1 eq $read2) && ($align1 eq $align2) && ($spec1 eq $spec2) && ($align1 ne "*") && ($spec1 eq $species1))
    {
      $geneHash{$gene1}{$spec1} += 1;
      open(OUT,">> $read_out") or die "Error writing to $read_out\n";
      print OUT "\n$read1\t$seq1\t$seq2\t$gene1\t$spec1\t$gene2\t$spec2\tGood";
      close OUT;
    }
    elsif(($read1 eq $read2) && ($align1 eq $align2) && ($spec1 eq $spec2) && ($align1 ne "*") && ($spec1 eq $species2))
    {
      if($sech2mel{$gene1})
      {
        $geneHash{$sech2mel{$gene1}}{$spec1} += 1;
        open(OUT,">> $read_out") or die "Error writing to $read_out\n";
        print OUT "\n$read1\t$seq1\t$seq2\t$gene1\t$spec1\t$gene2\t$spec2\tGood";
        close OUT;
      }
    }
    elsif(($read1 eq $read2) && ($align1 ne "*") && ($align2 eq "*") && ($spec1 eq $species1))
    {
      $geneHash{$gene1}{$spec1} += 1;
      open(OUT,">> $read_out") or die "Error writing to $read_out\n";
      print OUT "\n$read1\t$seq1\t$seq2\t$gene1\t$spec1\t$align2\t$align2\tOnly mate 1";
      close OUT;
    }
    elsif(($read1 eq $read2) && ($align1 ne "*") && ($align2 eq "*") && ($spec1 eq $species2))
    {
      if($sech2mel{$gene1})
      {
        $geneHash{$sech2mel{$gene1}}{$spec1} += 1;
        open(OUT,">> $read_out") or die "Error writing to $read_out\n";
        print OUT "\n$read1\t$seq1\t$seq2\t$gene1\t$spec1\t$align2\t$align2\tOnly mate 1";
        close OUT;
      }
    }
    elsif(($read1 eq $read2) && ($align1 eq "*") && ($align2 ne "*") && ($spec2 eq $species1))
    {
      $geneHash{$gene2}{$spec2} += 1;
      open(OUT,">> $read_out") or die "Error writing to $read_out\n";
      print OUT "\n$read1\t$seq1\t$seq2\t$align1\t$align1\t$gene2\t$spec2\tOnly mate 2";
      close OUT;
    }
    elsif(($read1 eq $read2) && ($align1 eq "*") && ($align2 ne "*") && ($spec2 eq $species2))
    {
      if($sech2mel{$gene2})
      {
        $geneHash{$sech2mel{$gene2}}{$spec2} += 1;
        open(OUT,">> $read_out") or die "Error writing to $read_out\n";
        print OUT "\n$read1\t$seq1\t$seq2\t$align1\t$align1\t$gene2\t$spec2\tOnly mate 2";
        close OUT;
      }
    }
    elsif(($read1 eq $read2) && ($align1 ne $align2) && ($spec1 ne $spec2) && ($spec1 eq $species1))
    {
      open(OUT,">> $read_out") or die "Error writing to $read_out\n";
      print OUT "\n$read1\t$seq1\t$seq2\t$gene1\t$spec1\t$gene2\t$spec2\tTransplicing";
      close OUT;
    }
    elsif(($read1 eq $read2) && ($align1 ne $align2) && ($spec1 ne $spec2) && ($spec1 eq $species2))
    {
      open(OUT,">> $read_out") or die "Error writing to $read_out\n";
      print OUT "\n$read1\t$seq1\t$seq2\t$gene1\t$spec1\t$gene2\t$spec2\tTransplicing";
      close OUT;
    }
    elsif(($read1 eq $read2) && ($align1 eq $align2) && ($align1 eq "*"))
    {
      open(OUT,">> $read_out") or die "Error writing to $read_out\n";
      print OUT "\n$read1\t$seq1\t$seq2\t$align1\t$align1\t$align2\t$align2\tNo alignment";
      close OUT;
    }
    else
    {
      open(OUT,">> $read_out") or die "Error writing to $read_out\n";
      print OUT "\n$read1\t$seq1\t$seq2\t$gene1\t$spec1\t$gene2\t$spec2\tInvestigate";
      close OUT;
    }
  }
}
close SAM1;
close SAM2;

my$gene_out = $ARGV[6];
open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
print OUT "gene\t$species1\t$species2";
close OUT;

foreach(sort keys %geneHash)
{
  if($geneHash{$_}{$species1} && $geneHash{$_}{$species2})
  {
    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}";
    close OUT;
  }
  elsif($geneHash{$_}{$species1})
  {
    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$_\t$geneHash{$_}{$species1}\t0";
    close OUT;
  }
  elsif($geneHash{$_}{$species2})
  {
    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$_\t0\t$geneHash{$_}{$species2}";
    close OUT;
  }
}
