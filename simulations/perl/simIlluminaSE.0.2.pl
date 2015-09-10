#!/usr/bin/perl

########################################################################
# 
# 05/13/2010
#
# simIlluminaSingle0.2.pl
#
# Purpose: simulate sets of single reads from reference databases
# 
# Input: a FASTA reference and SNP locations (with bp info)
#
# Output: millions of short sequence reads with various properties (mutation rates,
#         lengths, etc.)
#
# Fixed: new version 0.2 --> generate reads around SNPs (categorized as "both"), 
#                            Basically, generate reads from the entire reference
#
#        05/14/2010 --> created gene (or chrom) hash structure instead of going through file
#
# Notes: 05/13/2010 --> need to fix gene (or chromosome) loop so that we only
#                       go through the reference file once
#
########################################################################

use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

srand(12345);

my$length = $ARGV[0]; # read length
my$error = $ARGV[1]; # base-calling error rate

my$match; my$As; my$Cs; my$Gs; my$Ts; my$Ns;
my$line; my$gene; my$pos; my$cons; my$cov; my$call;
my@elements; my@alleles;
my%SNPhash;

# now create a hash of SNPs to be used later to generate reads around
my$SNP_file = $ARGV[2];
open(SNP,"$SNP_file") or die "\nError opening $SNP_file\n";
while(<SNP>)
{
  $match = 0;
  $As = 0;
  $Cs = 0;
  $Gs = 0;
  $Ts = 0; 
  $Ns = 0;

  chomp;  
  @elements = split(/\s+/,$_);
  $gene = $elements[0];
  $pos = $elements[1];
  $cons = $elements[3];
  $cov = $elements[7];
  $call = $elements[8];
  #print "Call:{$call}\t";
  @alleles = split("",$call);
  #foreach(@alleles){print "$_";}

  foreach(@alleles)
  {
    if(/[.,]/){$match+=1;}
    elsif(/[Aa]/){$As+=1;}
    elsif(/[Cc]/){$Cs+=1;}
    elsif(/[Gg]/){$Gs+=1;}
    elsif(/[Tt]/){$Ts+=1;}
    elsif(/[Nn]/){$Ns+=1;}
  }
  #print "$gene\t$pos\t$call\t$As\t$Cs\t$Gs\t$Ts\n";

  unless($cov<20 && (($As>0&&$Cs>0)||($As>0&&$Gs>0)||($As>0&&$Ts>0)||($Cs>0&&$Gs>0)||($Cs>0&&$Ts>0)||($Gs>0&&$Ts>0)||$Ns>0))
  {
    unless($match==$cov)
    {
      $SNPhash{$gene}{$pos} = $cons;
    }
  }    
}
close SNP;

my$seq; 
my%geneHash;

my$ref = $ARGV[3];
open(REF,"$ref") or die "\nError opening $ref\n";

while(<REF>)    
{
  chomp;
  if(/^>(\S+)/)
  {
    $gene = $1;
    undef $seq;
  }
  else
  {
    $seq = $seq . $_;
    $geneHash{$gene} = $seq;  
  }
}
close REF;

#foreach $gene (keys %geneHash)
#{
#  print "Gene: $gene\nSeq: $geneHash{$gene}\n";
#}

my$i; my$j; my$k = 1;
my$reflength; my$posNeg;
my$readRef; my$readCons; my$false; my$baseRef; my$baseCons;
my@nuc; my@nucNeg;

my$out = $ARGV[4];
open(OUT,">$out") or die "\nError opening $out\n";

foreach $gene (keys %geneHash) 
{
  @nuc = split("",$geneHash{$gene});
  @nucNeg = reverse @nuc;
  $reflength = scalar@nuc;
  #print "Length of reference: $reflength\n";

  foreach $pos (keys %{$SNPhash{$gene}})
  {
    unless($length > $reflength)
    {
      # positive orientation
      for($i=0;$i<$reflength-$length+1;$i++) # generate reads from the whole reference
      {
        for($j=0;$j<$length;$j++)
        {
          if($error == 0)
	  {
            $readRef = $readRef . $nuc[$i+$j];
	    if($i+$j == $pos-1){$readCons = $readCons . $SNPhash{$gene}{$pos};}
            else{$readCons = $readCons . $nuc[$i+$j];}
          }
          else
	  {
            $false = rand();
            #print "i = $i\tj = $j\n";
            if($false < $error)
            {
              $baseRef = rand();
              $baseCons = rand();
              #print "Error in base $j!\n";
              #print "Error: $false\tbaseRef: $baseRef\tbaseCons = $baseCons\n";
              if($baseRef < 0.25){$readRef = $readRef . "A";}
              elsif($baseRef >= 0.25 && $baseRef < 0.5){$readRef = $readRef . "C";}
              elsif($baseRef >= 0.5 && $baseRef < 0.75){$readRef = $readRef . "G";}
              else{$readRef = $readRef . "T";}

              if($baseCons < 0.25){$readCons = $readCons . "A";}
              elsif($baseCons >= 0.25 && $baseCons < 0.5){$readCons = $readCons . "C";}
              elsif($baseCons >= 0.5 && $baseCons < 0.75){$readCons = $readCons . "G";}
              else{$readCons = $readCons . "T";}
	    }
            else
	    {
              $readRef = $readRef . $nuc[$i+$j];
	      if($i+$j == $posNeg-1){$readCons = $readCons . $SNPhash{$gene}{$pos};}
              else{$readCons = $readCons . $nuc[$i+$j];}
            }
	  }
	}
        #print "i = $i, j = $j\n";
        if($pos-1 >= $i && $pos-1 <= $i+$length-1) # these are reads containing SNPs
	{
          print OUT "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_match\(\+\)\n$readRef\n";
          print OUT "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_mismatch\(\+\)\n$readCons\n";
        }
        else # these are reads NOT containing SNPs
	{
          print OUT "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\(\+\)\n$readRef\n";
          print OUT "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\(\+\)\n$readCons\n";
        }
        undef $readRef;
        undef $readCons;

        $k += 1;    
      }
      $k = 1;
     
      # negative orientation
      $posNeg = $reflength - $pos + 1; # reverse the position of the SNP
      for($i=0;$i<$reflength-$length+1;$i++) # generate reads across entire reference
      {
        for($j=0;$j<$length;$j++)
        {
          if($error == 0)
	  {
            $readRef = $readRef . $nucNeg[$i+$j];
	    if($i+$j == $posNeg-1){$readCons = $readCons . $SNPhash{$gene}{$pos};}
            else{$readCons = $readCons . $nucNeg[$i+$j];}
          }
          else
	  {
            $false = rand();
            #print "i = $i\tj = $j\n";
            if($false < $error)
            {
              $baseRef = rand();
              $baseCons = rand();
              #print "Error in base $j!\n";
              #print "Error: $false\tbaseRef: $baseRef\tbaseCons = $baseCons\n";
              if($baseRef < 0.25){$readRef = $readRef . "A";}
              elsif($baseRef >= 0.25 && $baseRef < 0.5){$readRef = $readRef . "C";}
              elsif($baseRef >= 0.5 && $baseRef < 0.75){$readRef = $readRef . "G";}
              else{$readRef = $readRef . "T";}

              if($baseCons < 0.25){$readCons = $readCons . "A";}
              elsif($baseCons >= 0.25 && $baseCons < 0.5){$readCons = $readCons . "C";}
              elsif($baseCons >= 0.5 && $baseCons < 0.75){$readCons = $readCons . "G";}
              else{$readCons = $readCons . "T";}
	    }
            else
	    {
              $readRef = $readRef . $nucNeg[$i+$j];
	      if($i+$j == $posNeg-1){$readCons = $readCons . $SNPhash{$gene}{$pos};}
              else{$readCons = $readCons . $nucNeg[$i+$j];}
            }
	  }
	}
        #print "i = $i, j = $j\n";          
        if($posNeg-1 >= $i && $posNeg-1 <= $i+$length-1) # these are reads containing SNPs
	{
          print OUT "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_match\(\-\)\n$readRef\n";
          print OUT "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_mismatch\(\-\)\n$readCons\n";
        }
        else # these are reads NOT containing SNPs
	{
          print OUT "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\(\-\)\n$readRef\n";
          print OUT "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\(\-\)\n$readCons\n";
        }
   
        undef $readRef;
        undef $readCons;

        $k += 1;
      }
      $k = 1;
    }  
  }
}
close OUT;



