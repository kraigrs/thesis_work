#!/usr/bin/perl

########################################################################
# 
# 10/04/2011
#
# edit_genome.pl
#
# Purpose: take in a reference genome and return the edited alternate genome and ambiguous genome
# 
# Input: SNP file, original reference, directory for new references 
#
# Output: outputs bad SNPs (SNP showing inconsistency with reference)
#         
# Syntax: perl edit_genome.pl ../Dmel/chromFa/chr2L.fa ../DGRP/chr2L_SNPs.txt ../Dmel/chromFa > ../DGRP/chr2L_badSNPs.txt 
#
########################################################################

use strict;
use warnings;
use Switch;

my$ref = $ARGV[0];
my$SNP_file = $ARGV[1];
my$ref_dir = $ARGV[2];

our($chr,$pos,$code,$seq,$i,$base,$j,$k);
our(@elements,@refPos);
our(%SNPhash);

open(SNP,"$SNP_file") or die "\nError opening $SNP_file\n";
while(<SNP>)
{
  chomp;  
  @elements = split(/\s+/,$_);
  $chr = $elements[0];
  $pos = $elements[1];
  $code = $elements[2];
  
  unless($code eq "N"){$SNPhash{$pos} = $code;}
    
}
close SNP;

#print "\ninputting reference\n";
#$i = 0;
open(REF,"$ref") or die "\nError opening $ref\n";
while(<REF>)    
{
  chomp;
  #print "Here!\n";
  if(/^>chr(\S+)$/){$chr = $1;}
  else
  {
    @elements = split("",$_); 
    foreach(@elements){push @refPos,uc($_);}
  }
  #$j = ++$i; print "$j\n"; 
}
close REF;
#print "\ndone\n\n";

#specify alternate genomes
open(ALT,">$ref_dir/chr$chr\_alt.fa") or die "\nError opening $ref_dir/chr$chr\_alt.fa\n";
open(AMB,">$ref_dir/chr$chr\_ambig.fa") or die "\nError opening $ref_dir/chr$chr\_ambig.fa\n";

print ALT "\>chr$chr\n";
print AMB "\>chr$chr\n";

$j = 0;
for($i=0;$i<@refPos;$i++) 
{
  if($refPos[$i] eq "N")
  {
    print ALT "N";
    print AMB "N";
    $j += 1;
  }
  elsif($SNPhash{$i+1})
  {
    $base = &altAllele($refPos[$i],$SNPhash{$i+1});
    if($base eq "N")
    {
      $k = $i+1; 
      print "$chr\t$k\t$SNPhash{$k}\n";
      print ALT "$refPos[$i]";
      print AMB "$refPos[$i]";
    }
    else
    {
      print ALT "$base";
      print AMB "$SNPhash{$i+1}";
    }
    $j += 1;
  }
  else
  {
    $base = $refPos[$i];
    print ALT "$base";
    print AMB "$base";
    $j += 1;  
  }
  
  if($j == 50){print ALT "\n"; print AMB "\n"; $j = 0;}
}
#print "j=$j\n";
unless($j == 0){print ALT "\n"; print AMB "\n";}

close ALT; close AMB;
#print "\nNumber mismatches between reference and DGRP SNP calls: $k\n\n";

sub altAllele
{
  my($n,$amb) = @_;
  my$var;
  switch ($amb)
  {
    case "R"
    {
      if($n eq "A"){$var = "G";}
      elsif($n eq "G"){$var = "A";}
      else{$var = "N";}
    }
    case "Y"
    {
      if($n eq "C"){$var = "T";}
      elsif($n eq "T"){$var = "C";}
      else{$var = "N";}
    }
    case "M"
    {
      if($n eq "A"){$var = "C";}
      elsif($n eq "C"){$var = "A";}
      else{$var = "N";}
    }
    case "K" 
    {
      if($n eq "G"){$var = "T";}
      elsif($n eq "T"){$var = "G";}
      else{$var = "N";}
    }
    case "S"
    {
      if($n eq "G"){$var = "C";}
      elsif($n eq "C"){$var = "G";}
      else{$var = "N";}
    }
    case "W"
    {
      if($n eq "A"){$var = "T";}
      elsif($n eq "T"){$var = "A";}
      else{$var = "N";}
    }
    else{die "Something is seriously wrong here!\n";}
  }
  return $var;
}
