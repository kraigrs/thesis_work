#!/usr/bin/perl

########################################################################
# 
# 11/12/2013
#
# vcf2fasta.pl
#
# Purpose: take in a VCF file and reference FASTA file and output the alternative reference
# 
# Input: VCF file, original reference, name of alternative reference 
#
# Output: alternative reference in FASTA format
#         
# Syntax: perl vcf2fasta.pl <VCF file> <FASTA reference> <name of alternative reference> 
#
########################################################################

use strict;
use warnings;
use List::Util qw(min max);

my$vcf = $ARGV[0];
my$fasta = $ARGV[1];
my$output = $ARGV[2];

my$pos; my$seq; my$i; my$base; my$first; my$len; my$ref_allele; my$alt_allele;
my$line; my$locus; my$ref; my$alt; my$flag;
my@elements; my@meta;
my%genome; my%variants;

# read in reference sequence
$seq=""; $first=1;
open(REF,"$fasta") or die "\nError opening $fasta\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)\s+/)
  {
    unless($first == 1){$genome{$locus} = $seq;} # the first time, don't save sequence in hash
    $locus = $1;
    $seq = "";
    $first = 0;
  }
  else
  {
    $seq = $seq.uc($_);
  }
}
$genome{$locus} = $seq; # this needs to be done since last entry will not be stored in hash
close REF;

#foreach $locus (keys %genome)
#{
#  $len = length($genome{$locus});
#  print "$locus\t$len\n";
#}
#exit;

open(VCF,"$vcf") or die "\nError opening $vcf\n";
while(<VCF>)
{
  chomp;
  $line = $_;
  unless(/^#/)
  {
    @elements = split("\t",$_);
    $locus = $elements[0];
    $pos = $elements[1];
    $ref = $elements[3];
    $alt = $elements[4];
    if($elements[7] =~ /\S*(AC=(\d+|\d+\,\d+));/){$flag = $2;}

    #print "$locus\t$pos\t$ref\t$alt\t$flag\n"; exit;
    $variants{$locus}{$pos} = [$ref,$alt,$flag];
  }
}
close VCF;

#print "FBgn0248127\t52\t";
#print "$variants{FBgn0248127}{52}[0]\t";
#print "$variants{FBgn0248127}{52}[1]\t";
#print "$variants{FBgn0248127}{52}[2]\n";
#exit;

#specify alternate genome

open(ALT,"> $output") or die "\nError opening $output\n";
foreach $locus (keys %genome)
{
  print ALT "\>$locus\n";
  $len = length($genome{$locus});

  for($i=0;$i<$len;$i++) 
  {
    $j = $i + 1;
    $ref_allele = substr $genome{$locus}, $i, 1;

    if($variants{$locus}{$j})
    {
      $ref = $variants{$locus}{$j}[0];
      $alt = $variants{$locus}{$j}[1];
      $flag = $variants{$locus}{$j}[2];

      if(length($ref)+length($alt) == 2) # SNVs
      {
        if($flag == 1){$code = ambig($ref,$alt);} # heterozygous
        elsif($flag == 2){$code = $variants{$locus}{$j}[1];} # homozygous
        else{die "Error in SNV definition!\n";}
        print ALT "$code";
      }

      elsif($alt =~ /\,/) # alleles not matching the reference
      {
        if(length($ref) == 1 && length($alt) == 3) # SNVs
	{
          @genos = split(/\,/,$alt);
          $code = ambig($genos[0],$genos[1]);
          print ALT "$code";
        }
      }



    }
    elsif($variants{$locus}{$i})
    {
      if($ref_allele ne $genos{$chr}{$i})
      { 
        print "$chr\t$i\t$one\t$ref_allele\t$genos{$chr}{$i}\n";
        print ALT "$genos{$chr}{$i}";
      }
      else
      {
        print ALT "$ref_allele";
      }
    }
    else
    {
      print ALT "$ref_allele";
    }
  }
}
close ALT;

sub ambig
{
  my($a1,$a2) = @_;
  my$var;
  if($a1.$a2 eq "AC" || $a1.$a2 eq "CA"){$var = "M";}
  elsif($a1.$a2 eq "AG" || $a1.$a2 eq "GA"){$var = "R";}
  elsif($a1.$a2 eq "AT" || $a1.$a2 eq "TA"){$var = "W";}
  elsif($a1.$a2 eq "CG" || $a1.$a2 eq "GC"){$var = "S";}
  elsif($a1.$a2 eq "CT" || $a1.$a2 eq "TC"){$var = "Y";}
  elsif($a1.$a2 eq "GT" || $a1.$a2 eq "TG"){$var = "K";}
  else{$var = "N";}
  return $var;
}
