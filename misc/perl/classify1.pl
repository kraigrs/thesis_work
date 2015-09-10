#!/usr/bin/perl

###################################################################################
#
# 02/12/2009
#
# classify1.pl
#
# Purpose: (adapted from readInfo.pl) classify reads (single and paired-end) as matching
#          or mismatching the reference genome
# 
# Input: two separate alignment files (mates 1 and mates 2) 
#
# Output: a list of genes and their respective match and mismatch counts, along with a 
#         summary of all the reads 
#
# E.g. perl classify1.pl ../mel_sec_data/Hyb.dmel-gene.SNPs.txt ../mel_sec_data/Hyb_mate1.dmel-gene.sam ../mel_sec_data/Hyb_mate2.dmel-gene.sam ../mel_sec_data/Hyb.dmel-gene.reads.txt ../mel_sec_data/Hyb.dmel-gene.genes.txt &
#
###################################################################################

use strict;
use warnings;

my$match; my$As; my$Cs; my$Gs; my$Ts; my$Ns;
my$line; my$gene; my$pos; my$ref; my$cov; my$call;
my@elements; my@alleles;
my%SNPhash;

my$SNP_file = $ARGV[0];
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
  $line = $_;
  @elements = split(/\s/,$line);
  $gene = $elements[0];
  $pos = $elements[1];
  $ref = $elements[2];
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

  unless($cov<5 && (($As>0&&$Cs>0)||($As>0&&$Gs>0)||($As>0&&$Ts>0)||($Cs>0&&$Gs>0)||($Cs>0&&$Ts>0)||($Gs>0&&$Ts>0)||$Ns>0))
  {
    unless($match==$cov)
    {
      $SNPhash{$gene}{$pos} = $ref;
    }
  }    
}
close SNP;
print "\nFinished making SNP hash!\n";

my$line1; my$line2;
my$read1; my$read2;
my$gene1; my$gene2;
my$seq1; my$seq2;
my$m1start; my$m2start;
my@elements1; my@elements2;
my@seq1; my@seq2;
my$posKey;
my%geneHash;

my$mel1Ct = 0; my$sim1Ct = 0;
my$mel2Ct = 0; my$sim2Ct = 0;

my$m = 0; my$s = 0; my$b = 0;
my$flag;
my$MEL = 0; my$SIM = 0;

my$SAMfile1 = $ARGV[1];
my$SAMfile2 = $ARGV[2];

my$read_out = $ARGV[3]; 

open (SAM1, "$SAMfile1") or die "Can't open $SAMfile1 for reading!";
open (SAM2, "$SAMfile2") or die "Can't open $SAMfile2 for reading!";

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
    $gene1 = $elements1[2];
    $m1start = $elements1[3]; 
    $seq1 = $elements1[9];
  
    @elements2 = split(/\s/,$line2);
    if($elements2[0] =~ /^(HWI.+)\/\d{1}/){$read2 = $1;}
    $gene2 = $elements2[2];
    $m2start = $elements2[3]; 
    $seq2 = $elements2[9];

    if(($read1 eq $read2) && ($gene1 eq $gene2) && ($gene1 ne "*"))
    {
      $mel1Ct = 0;
      $sim1Ct = 0;
      $mel2Ct = 0;
      $sim2Ct = 0;
      $m = 0;
      $s = 0;
      $b = 0;
      $flag = 'no-call';
      $MEL = 0;
      $SIM = 0;
      @seq1 = split("",$seq1);
      #print "@seq1\n";
      @seq2 = split("",$seq2);
      #print "@seq2\n";
      for $posKey (keys %{$SNPhash{$gene1}}) 
      {
        #print "$posKey\n";
        if($posKey >= $m1start && $posKey <= ($m1start+scalar@seq1-1))
        {
          if($seq1[$posKey-$m1start] eq $SNPhash{$gene1}{$posKey}){$mel1Ct += 1;}
          else{$sim1Ct += 1;}
        }
        elsif($posKey >= $m2start && $posKey <= ($m2start+scalar@seq2-1))
        {
          if($seq2[$posKey-$m2start] eq $SNPhash{$gene1}{$posKey}){$mel2Ct += 1;}
          else{$sim2Ct += 1;} 
        } 
      }

      $m = $mel1Ct + $mel2Ct;   
      $s = $sim1Ct + $sim2Ct;  

      if($m + $s == 0)
      {
        $geneHash{$gene1}{'m'} += 0; 
        $geneHash{$gene1}{'s'} += 0;  
        $flag = 'no_SNPs';
        $b = 1;
      }
      elsif($m > 0 && $s == 0)
      {
        $geneHash{$gene1}{'m'} += 1; 
        $geneHash{$gene1}{'s'} += 0;
        $flag = 'confident';
        $MEL = 1;
      }
      elsif($m == 0 && $s > 0)
      {
        $geneHash{$gene1}{'m'} += 0; 
        $geneHash{$gene1}{'s'} += 1;     
        $flag = 'confident';
        $SIM = 1;
      }
      elsif(($mel1Ct > 1 && $sim1Ct == 1) && ($mel2Ct == 0 && $sim2Ct == 0))
      {
        $geneHash{$gene1}{'m'} += 0; 
        $geneHash{$gene1}{'s'} += 0;  
        $flag = 'sequencing_error';
      }
      elsif(($mel1Ct == 0 && $sim1Ct == 0) && ($mel2Ct > 1 && $sim2Ct == 1))
      {
        $geneHash{$gene1}{'m'} += 0; 
        $geneHash{$gene1}{'s'} += 0;  
        $flag = 'sequencing_error';
      } 
      elsif(($mel1Ct == 1 && $sim1Ct > 1) && ($mel2Ct == 0 && $sim2Ct == 0))
      {
        $geneHash{$gene1}{'m'} += 0; 
        $geneHash{$gene1}{'s'} += 0;  
        $flag = 'sequencing_error';
      }
      elsif(($mel1Ct == 0 && $sim1Ct == 0) && ($mel2Ct == 1 && $sim2Ct > 1))
      {
        $geneHash{$gene1}{'m'} += 0; 
        $geneHash{$gene1}{'s'} += 0;  
        $flag = 'sequencing_error';
      }
      elsif(($mel1Ct >= 1 && $sim1Ct == 0) && ($mel2Ct == 0 && $sim2Ct >= 1))
      {
        $geneHash{$gene1}{'m'} += 1; 
        $geneHash{$gene1}{'s'} += 1; 
        $flag = 'trans-splicing';
        $MEL = 1;
        $SIM = 1;
      }
      elsif(($mel1Ct == 0  && $sim1Ct >= 1) && ($mel2Ct >= 1 && $sim2Ct == 0))
      {
        $geneHash{$gene1}{'m'} += 1;
        $geneHash{$gene1}{'s'} += 1; 
        $flag = 'trans-splicing';
        $MEL = 1;
        $SIM = 1;
      }
      elsif( (($mel1Ct > 0 && $sim1Ct < 2) && ($mel2Ct > 0 && $sim2Ct == 0)) || (($mel1Ct > 0 && $sim1Ct == 0) && ($mel2Ct > 0 && $sim2Ct < 2)) )
      {
        $geneHash{$gene1}{'m'} += 0; 
        $geneHash{$gene1}{'s'} += 0;
        $flag = 'sequencing_error';
      }
      elsif( (($mel1Ct < 2 && $sim1Ct > 0) && ($mel2Ct == 0 && $sim2Ct > 0)) || (($mel1Ct == 0 && $sim1Ct > 0) && ($mel2Ct <= 1 && $sim2Ct > 0)) )
      {
        $geneHash{$gene1}{'m'} += 0; 
        $geneHash{$gene1}{'s'} += 0;
        $flag = 'sequencing_error';
      }
      else
      {
        $geneHash{$gene1}{'m'} += 0;
        $geneHash{$gene1}{'s'} += 0;
        $flag = 'no_call';
      }
      $geneHash{$gene1}{'b'} += $b;

      open(OUT,">> $read_out") or die "Error writing to $read_out\n";
      print OUT "\n$read1\t$gene1\t$seq1\t$seq2\t$mel1Ct\t$sim1Ct\t$mel2Ct\t$sim2Ct\t$b\t$MEL\t$SIM\t$flag";
      close OUT;
    }
    elsif(($read1 eq $read2) && ($gene1 ne "*") && ($gene2 eq "*"))
    {
      $mel1Ct = 0;
      $sim1Ct = 0;
      $m = 0;
      $s = 0;
      $b = 0;
      $flag = 'no-call';
      $MEL = 0;
      $SIM = 0;
      @seq1 = split("",$seq1);
      #print "@seq1\n";
      for $posKey (keys %{$SNPhash{$gene1}}) 
      {
        if($posKey >= $m1start && $posKey <= ($m1start+scalar@seq1-1))
        {
          if($seq1[$posKey-$m1start] eq $SNPhash{$gene1}{$posKey}){$mel1Ct += 1;}
          else{$sim1Ct += 1;}
        }
      }
      if($mel1Ct + $sim1Ct == 0)
      {
        $geneHash{$gene1}{'m'} += 0; 
        $geneHash{$gene1}{'s'} += 0;  
        $flag = 'no_SNPs';
        $b = 1;
      }
      elsif($mel1Ct > 0 && $sim1Ct == 0)
      {
        $geneHash{$gene1}{'m'} += 1; 
        $geneHash{$gene1}{'s'} += 0;
        $flag = 'confident';
        $MEL = 1;
      }
      elsif($mel1Ct == 0 && $sim1Ct > 0)
      {
        $geneHash{$gene1}{'m'} += 0; 
        $geneHash{$gene1}{'s'} += 1;     
        $flag = 'confident';
        $SIM = 1;
      }
      open(OUT,">> $read_out") or die "Error writing to $read_out\n";
      print OUT "\n$read1\t$gene1\t$seq1\tNA\t$mel1Ct\t$sim1Ct\tNA\tNA\t$b\t$MEL\t$SIM\t$flag";
      close OUT;
    }
    elsif(($read1 eq $read2) && ($gene1 eq "*") && ($gene2 ne "*"))
    {
      $mel2Ct = 0;
      $sim2Ct = 0;
      $m = 0;
      $s = 0;
      $b = 0;
      $flag = 'no-call';
      $MEL = 0;
      $SIM = 0;
      for $posKey (keys %{$SNPhash{$gene2}}) 
      {
        #print "$posKey\n";
        if($posKey >= $m2start && $posKey <= ($m2start+scalar@seq2-1))
        {
          if($seq2[$posKey-$m2start] eq $SNPhash{$gene2}{$posKey}){$mel2Ct += 1;}
          else{$sim2Ct += 1;} 
        } 
      }
      if($mel2Ct + $sim2Ct == 0)
      {
        $geneHash{$gene2}{'m'} += 0; 
        $geneHash{$gene2}{'s'} += 0;  
        $flag = 'no_SNPs';
        $b = 1;
      }
      elsif($mel2Ct > 0 && $sim2Ct == 0)
      {
        $geneHash{$gene2}{'m'} += 1; 
        $geneHash{$gene2}{'s'} += 0;
        $flag = 'confident';
        $MEL = 1;
      }
      elsif($mel2Ct == 0 && $sim2Ct > 0)
      {
        $geneHash{$gene2}{'m'} += 0; 
        $geneHash{$gene2}{'s'} += 1;     
        $flag = 'confident';
        $SIM = 1;
      }
      open(OUT,">> $read_out") or die "Error writing to $read_out\n";
      print OUT "\n$read2\t$gene2\tNA\t$seq2\tNA\tNA\t$mel2Ct\t$sim2Ct\t$b\t$MEL\t$SIM\t$flag";
      close OUT;
    }
  }
}
close SAM1;
close SAM2;

my$madj;
my$sadj;

my$gene_out = $ARGV[4];
open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
print OUT "gene\tmel\tsim\tboth\tmelAdj\tsimAdj";
close OUT;

for my$k1 (sort keys %geneHash)
{
  if(($geneHash{$k1}{'m'} > 0 || $geneHash{$k1}{'s'} > 0) && $k1 ne "*")
  {
    $madj = $geneHash{$k1}{'m'} + ($geneHash{$k1}{'m'} / ($geneHash{$k1}{'m'} + $geneHash{$k1}{'s'})) * $geneHash{$k1}{'b'};  
    $sadj = $geneHash{$k1}{'s'} + ($geneHash{$k1}{'s'} / ($geneHash{$k1}{'m'} + $geneHash{$k1}{'s'})) * $geneHash{$k1}{'b'};
    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$k1\t$geneHash{$k1}{m}\t$geneHash{$k1}{s}\t$geneHash{$k1}{b}\t$madj\t$sadj";
    close OUT;
  }
}
