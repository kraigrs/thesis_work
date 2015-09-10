#!/usr/bin/perl

########################################################################
# 
# 03/22/2010
#
# simIlluminaSingle0.1.pl
#
# Purpose: simulate sets of single reads from reference databases
# 
# Input: a FASTA reference and SNP locations (with bp info)
#
# Output: millions of short sequence reads with various properties (mutation rates,
#         lengths, etc.)
#
# Fixed: 05/11/2010 --> replaced last 3 loops with a single looping structure, more efficient 
#                   --> should consider replacing "find gene" loop with making a gene hash, 
#                       but may require too much memory 
#
#        05/12/2010 --> added both strand orientation simulated reads (basically doubling coverage)
#
# Notes: 05/12/2010 --> may be beneficial to also simulate reads other than those around SNPs.
#                       This may better represent reads that we later classify as "both".
#                       The simulation strategy would be something like this:
#                       1). Choose a gene (or chromosome) to start with, and generate reads starting 
#                           at the beginning of the reference sequence
#                       2). When a SNP is reached, start original process.
#
#        05/12/2010 --> also need to simulate reads in the opposite orientation!!!
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

my$ref = $ARGV[3];
my$out = $ARGV[4];

open(OUT,">$out") or die "\nError opening $out\n";

my$seq; my$found = 0; my$i; my$j; my$k = 1;
my$reflength; my$pospluslength; my$posNeg;
my$readRef; my$readCons; my$false; my$baseRef; my$baseCons;
my@nuc; my@nucNeg;

foreach $gene (sort keys %SNPhash)
{
  open(REF,"$ref") or die "\nError opening $ref\n";
  HERE: while(<REF>)    
  {
    chomp;
    if(/^>$gene\s+/)
    {      
      $found = 1; 
    }
    elsif(/^>.+/ && $found)
    {   
      last HERE;
    }
    elsif($found)
    {
      $seq = $seq . $_;  
    }
  }
  close REF;
  $found = 0;
  #print "$gene\n$seq\n";
  @nuc = split("",$seq);
  @nucNeg = reverse @nuc;
  undef $seq;
  $reflength = scalar@nuc;
  #print "Length of reference: $reflength\n";

  foreach $pos (sort keys %{$SNPhash{$gene}})
  {
    unless($length > $reflength)
    {
      # positive orientation
      for($i=max($pos-$length,0);$i<min($reflength-$length+1,$pos);$i++) # account for read exceeding reference
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
        print OUT "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_match\(\+\)\n$readRef\n";
        print OUT "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_mismatch\(\+\)\n$readCons\n";

        undef $readRef;
        undef $readCons;

        $k += 1;    
      }
      $k = 1;
     
      # negative orientation
      $posNeg = $reflength - $pos + 1; # reverse the position of the SNP
      for($i=max($posNeg-$length,0);$i<min($reflength-$length+1,$posNeg);$i++) # account for read exceeding reference
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
          
        print OUT "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_match\(\-\)\n$readRef\n";
        print OUT "\>HWI_Gene_$gene\_SNP_$pos\_Read_$k\_mismatch\(\-\)\n$readCons\n";
   
        undef $readRef;
        undef $readCons;

        $k += 1;
      }
      $k = 1;
    }  
  }
}
close OUT;

# alternative gene hash structure

#open(REF,"$ref") or die "\nError opening $ref\n";
#while(<REF>)    
#{
#  chomp;
#  if(/^>(\S+)/)
#  {
#    $gene = $1;
#    undef $seq;
#  }
#  else
#  {
#    $seq = $seq . $_;
#    $geneHash{$gene} = $seq;  
#  }
#}
#close REF;


