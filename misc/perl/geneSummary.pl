#!/usr/bin/perl

###################################################################################
#
# 11/12/2009
#
# geneSummary.pl
#
# Purpose: from the SNP output from SAMtools, count SNPs per gene. Also, count
#          # of reads on a gene from .sam output
# 
# Input: text file containing all the SNPs within a gene and .sam file with all
#        aligned reads and genes
#
# Output: list of genes and their SNP counts and read counts
#
###################################################################################

use strict;
use warnings;

my$starttime = time;

my$SNP_file = $ARGV[0];

my$match;
my$As;
my$Cs;
my$Gs;
my$Ts;
my$Ns;
my$mismatch;
my@elements;
my$gene;
my$call;
my@alleles;
my%SNPhash;

open(SNP,"$SNP_file") or die "\nError opening $SNP_file\n";
while(<SNP>)
{
  $match = 0;
  $As = 0;
  $Cs = 0;
  $Gs = 0;
  $Ts = 0; 
  $Ns = 0;
  $mismatch = 0;

  chomp $_;  
 
  @elements = split(/\s/,$_);
  $gene = $elements[0];
  $call = $elements[8];
  #print "Call:{$call}\t";
  @alleles = split("",$call);
  #foreach(@alleles){print "$_";}
  foreach(@alleles)
  {
    if(/[\.,]/){$match+=1;}
    elsif(/[Aa]/){$As+=1;}
    elsif(/[Cc]/){$Cs+=1;}
    elsif(/[Gg]/){$Gs+=1;}
    elsif(/[Tt]/){$Ts+=1;}
    elsif(/[Nn]/){$Ns+=1;}
  }
  $mismatch = $As+$Cs+$Gs+$Ts;
  #print "$gene\t$pos\t$call\t$As\t$Cs\t$Gs\t$Ts\n";
  unless((($As>0&&$Cs>0)||($As>0&&$Gs>0)||($As>0&&$Ts>0)||($Cs>0&&$Gs>0)||($Cs>0&&$Ts>0)||($Gs>0&&$Ts>0))||$Ns>0)
  {
    $SNPhash{$gene}{'match'} += $match;
    $SNPhash{$gene}{'mismatch'} += $mismatch;
    $SNPhash{$gene}{'SNPs'} += 1;
  }    
}
close SNP;
print "\nFinished making SNP hash!\n";

my$SAM_file = $ARGV[1];

my$mate1found;
my$mate2found;
my$read1;
my$read2;
my$gene1;
my$gene2;
my%SAMhash;

print "\nParsing SAM output...\n";
open(SAM,"$SAM_file") or die "\nError opening $SAM_file\n";
while(<SAM>)
{
  chomp $_;  
  if(/^HWI.+\/1.+/)
  {
    $mate1found = 1;
    @elements = split(/\s/,$_);
    #print "$elements[0]";
    #print "$_\n";
    if($elements[0] =~ /^(HWI.+)\/\d{1}/){$read1 = $1;}
    #print "mate 1 found: $read1\n";
    $gene1 = $elements[2];   
  }
  
  elsif(/^HWI.+\/2.+/)
  {
    $mate2found = 1;
    @elements = split(/\s/,$_);
    #print "$elements[0]";
    #print "$_\n";
    if($elements[0] =~ /^(HWI.+)\/\d{1}/){$read2 = $1;}
    #print "mate 2 found: $read2\n";
    $gene2 = $elements[2];
  }
  
  if($mate1found && $mate2found && ($read1 eq $read2) && ($gene1 eq $gene2)){$SAMhash{$gene1} += 1;}
}

#for my$k1 (sort keys %SNPhash)
#{
#  print "$k1\t$SNPhash{$k1}\t$SAMhash{$k1}\n"; 
#}

my$out_file = $ARGV[2];

open(OUT,">> $out_file") or die "Error writing to $out_file\n";
print OUT "gene\tSNPs\treads\tmatches\tmismatches";
close OUT;

foreach(keys %SAMhash)
{
  if($SNPhash{$_})
  {
    open(OUT,">> $out_file") or die "Error writing to $out_file\n";
    print OUT "\n$_\t$SNPhash{$_}{SNPs}\t$SAMhash{$_}\t$SNPhash{$_}{match}\t$SNPhash{$_}{mismatch}";
    close OUT;
  }
  else
  {
    open(OUT,">> $out_file") or die "Error writing to $out_file\n";
    print OUT "\n$_\t0\t$SAMhash{$_}\t0\t0";
    close OUT;
  }
}

print "Done!\n";
printf ("\nTime elapsed: %d\n\n",time-$starttime);
