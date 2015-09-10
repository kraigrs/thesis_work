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

my$fasta = $ARGV[0];
my$vcf = $ARGV[1];
my$output = $ARGV[2];

my$locus; my$seq; my$i; my$pos; my$first; my$len; my$ref_allele;
my$line; my$ref; my$alt; my$flag; my$j; my$code;
my@elements; my@meta; my@genos;
my%genome; my%variants; my%genos;

# read in reference sequence
$seq=""; $first=1;
open(REF,"$fasta") or die "\nError opening $fasta\n";
while(<REF>)    
{
  chomp;
  if(/^\>(\S+)\s*/)
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
    @elements = split(/\s+/,$_);
    $locus = $elements[0];
    $pos = $elements[1];
    $ref = $elements[3];
    $alt = $elements[4];

    unless($elements[6] =~ /FAIL/)
    {      
      if($elements[7] =~ /\S*(AF=(1\.00|0\.500|0\.500,0\.500));/)
      {
        #print "$2\n";
        if($2 eq "1\.00")
        {
          $flag = "homo";
          #print "$line\n";
        }
        else
        {
          $flag = "hetero";
          #print "$line\n";
        }
      }

      #print "$locus\t$pos\t$ref\t$alt\t$flag\n"; exit;
      $variants{$locus}{$pos} = [$ref,$alt,$flag];
    }
  }
}
close VCF;

#$j = 0;
#foreach $locus (keys %variants)
#{
#  foreach $pos (keys %{$variants{$locus}})
#  {
#    if($variants{$locus}{$pos}[2] eq "hetero"){$j += 1;}
#  }
#}
#print "Total het variants being considered: $j\n"; exit;

# find variants

open(ALT,"> $output") or die "\nError opening $output\n";
foreach $locus (keys %genome)
{
  print ALT "\>$locus\n";
  $len = length($genome{$locus});
  #print "$locus\n";

  for($i=0;$i<$len;$i++) 
  {
    $j = $i+1;
    #print "$j\n";
    $ref_allele = substr $genome{$locus}, $i, 1;

    if($variants{$locus}{$j})
    {
      $ref = $variants{$locus}{$j}[0];
      $alt = $variants{$locus}{$j}[1];
      $flag = $variants{$locus}{$j}[2];

      if($alt =~ /\S+\,\S+/) # alleles not matching the reference
      {
        if(length($alt) == 3) # SNVs
	{
          #@genos = split(/\,/,$alt);
          #$code = ambig($genos[0],$genos[1]);
          $code = "N";
          print ALT "$code";
        }
        else
	{
          #$code = het_indels($ref,$alt);
          #if(length($code) == 1){print ALT "$code";}
          if(length($ref) > 1)
	  {
            $code = "N" x length($ref);
            $i += length($ref) - 1; # jump ahead since replacing bases in the ref, treat like deletion
            print ALT "$code";
          }
          else
          { 
            $code = "N";
            print ALT "$code";
	  }
	}
      }
      elsif(length($ref)+length($alt) == 2) # SNVs
      {
        if($ref ne $ref_allele){print "Error: $locus\t$ref_allele\t$ref\t$alt\n";}

        if($flag eq "hetero")
        {
          #$code = ambig($ref,$alt);
          $code = "N";
        }
        elsif($flag eq "homo"){$code = $alt;}
        else{die "Error in SNV definition!\n";}

        print ALT "$code";
      }
      else # indels
      {
        if($flag eq "hetero")
        {
          #$code = het_indels($ref,$alt);
          if(length($ref) > 1)
	  {
            $code = "N" x length($ref);
            $i += length($ref) - 1; # jump ahead since replacing bases in the ref, treat like deletion
            print ALT "$code";
          }
          else
          { 
            $code = "N";
            print ALT "$code";
	  }
        } 
        elsif($flag eq "homo")
        {
          $code = $alt;
          if(length($ref) > 1 && length($ref) > length($alt)) # deletion
	  {
            $i += length($ref)-1; # move ahead bases for deletion
            print ALT "$code";
	  }
          elsif(length($ref) < length($alt)) # insertion
          {
            if(length($ref) == 1){print ALT "$code";}
            else
	    {
              $i += length($ref)-1; # move ahead bases
              print ALT "$code";
	    }
	  }
        }
        else{die "Error in indel definition!\n";}
      }
    }
    else{print ALT "$ref_allele";}
  }
  print ALT "\n";
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

sub het_indels
{
  my($a1,$a2) = @_;
  my$n; my$m;
  my@vars;
  my@sizes;
  if($a2 =~ /\S+\,\S+/){@vars = split(/\,/,$a2);}
  else{$vars[0] = $a2;}
  push @vars,$a1;
  foreach(@vars){push @sizes, length($_);}
  $n = max @sizes;
  $m = "N" x $n;
  return $m;
}
