#!/usr/bin/perl

###################################################################################
#
# 12/01/2010
#
# classify.0.4.pl
#
# Purpose: combine 4 .bed files (lifted) to elucidate allele-specific expression 
# 
# Input: two separate alignment files (1 for each species)
#
# Output: a list of genes and their respective match and mismatch counts, along with a 
#         summary of all the reads 
#
# Fixes: this script represents an updated version fixing a memory leak problem. Also, instead of repeatedly
#        building gene and exon hashes, those are made once in the beginning
#
# E.g. perl classify.0.4.pl <const> <species1> <species2> <bedmate1species1> <bedmate2species1> <bedmate1species2> <bedmate2species2> <readout> <exonout> <geneout>
#     
#      <const>            ==> list of constitutive exons
#      <species1>         ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <species2>         ==> indicate name of 2nd species aligned to
#      <bedmate1species1> ==> file with mate 1 alignments to 1st species
#      <bedmate2species1> ==> file with mate 2 alignments to 1st species
#      <bedmate1species2> ==> file with mate 1 alignments to 2nd species
#      <bedmate2species2> ==> file with mate 2 alignments to 2nd species
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
my$line;
my$gene;
my$exon;
my$read;
my$name;
my@elements;
my%geneHash;
my%exonHash;
my%readList;

my$const = $ARGV[0];
my$species1 = $ARGV[1]; # Dmel
my$species2 = $ARGV[2]; # Dsec
my$bedmate1species1 = $ARGV[3];
my$bedmate2species1 = $ARGV[4];
my$bedmate1species2 = $ARGV[5];
my$bedmate2species2 = $ARGV[6];
my$out_dir = $ARGV[7];
my$suffix = $ARGV[8];

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
    $exon = $elements[1] . "," . $elements[2];

    $geneHash{$gene}{$species1} = 0; 
    $exonHash{$gene}{$exon}{$species1} = 0;

    $geneHash{$gene}{$species2} = 0; 
    $exonHash{$gene}{$exon}{$species2} = 0;

    $geneHash{$gene}{b} = 0; 
    $exonHash{$gene}{$exon}{b} = 0;
  }
}
close CONST;
print "\nDone with genes and exons hashes!\n";

# mate 1, species 1
open(BED,"$bedmate1species1") or die "Can't open $bedmate1species1 for reading!";
LINE:while(<BED>) 
{
  chomp;
  $line = $_;  
  #print "new line\n";
  if($line =~ /^\@.+/){next LINE;}
  elsif($line =~ /HWI.+\/\d+/)
  { 
    @elements = split(/\s/,$line);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene = $elements[0];
      $exon = $elements[1] . "," . $elements[2];
      $name = "1_1_".$read;
      #print "$name\t$gene\t$exon\n";

      if($readList{$name})
      {
        if($readList{$name}[0] eq $gene && $readList{$name}[1] eq $exon){next LINE;}
        else{delete $readList{$name}; next LINE;}
      }
      $readList{$name} = [$gene,$exon];
      #$readList{$name} = $gene;

      #print "$name\t$readList{$name}\n";
      #print "$name\t$readList{$name}[0]\t$readList{$name}[1]\n";
      #if($readList{'1_1_'.$read}){print "oh yeah!!!\n";}
      #print "$readList{'1_1_'.$read}[0]\t$readList{'1_1_'.$read}[1]\n";
    }
  }
} 
close BED;

print "\nDone with 1st hash!\n";

# mate 2, species 1
open(BED,"$bedmate2species1") or die "Can't open $bedmate2species1 for reading!";
LINE:while(<BED>) 
{
  chomp;
  $line = $_;  
  #print "new line\n";
  if($line =~ /^\@.+/){next LINE;}
  elsif($line =~ /HWI.+\/\d+/)
  { 
    @elements = split(/\s/,$line);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene = $elements[0];
      $exon = $elements[1] . "," . $elements[2];
      $name = "2_1_".$read;

      if($readList{$name})
      {
        if($readList{$name}[0] eq $gene && $readList{$name}[1] eq $exon){next LINE;}
        else{delete $readList{$name}; next LINE;}
      }
      $readList{$name} = [$gene,$exon];
    }
  }
}
close BED;

print "\nDone with 2nd hash!\n";

# mate 1, species 2
open(BED,"$bedmate1species2") or die "Can't open $bedmate1species2 for reading!";
LINE:while(<BED>) 
{
  chomp;
  $line = $_;  
  #print "new line\n";
  if($line =~ /^\@.+/){next LINE;}
  elsif($line =~ /HWI.+\/\d+/)
  { 
    @elements = split(/\s/,$line);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene = $elements[0];
      $exon = $elements[1] . "," . $elements[2];
      $name = "1_2_".$read;

      if($readList{$name})
      {
        if($readList{$name}[0] eq $gene && $readList{$name}[1] eq $exon){next LINE;}
        else{delete $readList{$name}; next LINE;}
      }
      $readList{$name} = [$gene,$exon];
    }
  }
}
close BED;

print "\nDone with 3rd hash!\n";

# mate 2, species 2
open(BED,"$bedmate2species2") or die "Can't open $bedmate2species2 for reading!";
LINE:while(<BED>) 
{
  chomp;
  $line = $_;  
  #print "new line\n";
  if($line =~ /^\@.+/){next LINE;}
  elsif($line =~ /HWI.+\/\d+/)
  { 
    @elements = split(/\s/,$line);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene = $elements[0];
      $exon = $elements[1] . "," . $elements[2];
      $name = "2_2_".$read;

      if($readList{$name})
      {
        if($readList{$name}[0] eq $gene && $readList{$name}[1] eq $exon){next LINE;}
        else{delete $readList{$name}; next LINE;}
      }
      $readList{$name} = [$gene,$exon];
    }
  }
}
close BED;

print "\nDone with 4th and final hash!\n";

##### start counting #####
#there should be 16 possibilities for alignments of both mates to both species (4 files, 2 possibilities each, 2^4 = 16)
#these can be taken care of on a case by case basis

open(OUT1,">$out_dir\/$suffix\.mosaik.reads.txt") or die "Error writing to $out_dir\/$suffix\.mosaik.reads.txt\n";
print OUT1 "read\tgene1$species1\tgene1$species2\tgene2$species1\tgene2$species2\tcall"; 

foreach $name (keys %readList)   
{
  @elements = split("_",$name);
  $read = $elements[2];   

  if($readList{'1_1_'.$read} && $readList{'1_2_'.$read}) # mate1 aligned to both genomes
  {
    if($readList{'2_1_'.$read} &&          #11
       $readList{'2_2_'.$read} && 
       $readList{'1_1_'.$read}[0] eq $readList{'2_1_'.$read}[0] && 
       $readList{'1_2_'.$read}[0] eq $readList{'2_2_'.$read}[0]) 
    {
      if($readList{'1_1_'.$read}[0] eq $readList{'1_2_'.$read}[0])
      {
        $geneHash{$readList{'1_1_'.$read}[0]}{b} += 1; 
        if($readList{'1_1_'.$read}[1] eq $readList{'2_1_'.$read}[1])
        {
          $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'1_1_'.$read}[1]}{b} += 1;
        } 
        else
        {
          $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'1_1_'.$read}[1]}{b} += 1; 
          $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'2_1_'.$read}[1]}{b} += 1;
        }     
        print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t$readList{'1_2_'.$read}[0]\t$readList{'2_1_'.$read}[0]\t$readList{'2_2_'.$read}[0]\tBoth";
      } 
      else{print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t$readList{'1_2_'.$read}[0]\t$readList{'2_1_'.$read}[0]\t$readList{'2_2_'.$read}[0]\tError11";}
    } 
    elsif($readList{'2_1_'.$read}) #1
    {
      if($readList{'1_1_'.$read}[0] eq $readList{'2_1_'.$read}[0])
      {
        $geneHash{$readList{'1_1_'.$read}[0]}{$species1} += 1; 
        if($readList{'1_1_'.$read}[1] eq $readList{'2_1_'.$read}[1])
        {
          $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'1_1_'.$read}[1]}{$species1} += 1;
        } 
        else
        {
          $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'1_1_'.$read}[1]}{$species1} += 1; 
          $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'2_1_'.$read}[1]}{$species1} += 1;
        } 
        print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t$readList{'1_2_'.$read}[0]\t$readList{'2_1_'.$read}[0]\t\*\t$species1";
      }
      else{print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t$readList{'1_2_'.$read}[0]\t$readList{'2_1_'.$read}[0]\t\*\tError1";}
    }
    elsif($readList{'2_2_'.$read}) #6
    {
      if($readList{'1_2_'.$read}[0] eq $readList{'2_2_'.$read}[0])
      {
        $geneHash{$readList{'1_2_'.$read}[0]}{$species2} += 1; 
        if($readList{'1_2_'.$read}[1] eq $readList{'2_2_'.$read}[1])
        {
          $exonHash{$readList{'1_2_'.$read}[0]}{$readList{'1_2_'.$read}[1]}{$species2} += 1;
        } 
        else 
        {
          $exonHash{$readList{'1_2_'.$read}[0]}{$readList{'1_2_'.$read}[1]}{$species2} += 1; 
          $exonHash{$readList{'1_2_'.$read}[0]}{$readList{'2_2_'.$read}[1]}{$species2} += 1;
        } 
        print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t$readList{'1_2_'.$read}[0]\t\*\t$readList{'2_2_'.$read}[0]\t$species2";
      }
      else{print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t$readList{'1_2_'.$read}[0]\t\*\t$readList{'2_2_'.$read}[0]\tError6";}
    }
    else #13
    {
      if($readList{'1_1_'.$read}[0] eq $readList{'1_2_'.$read}[0] && 
         $readList{'1_1_'.$read}[1] eq $readList{'1_2_'.$read}[1])
      {
        $geneHash{$readList{'1_1_'.$read}[0]}{b} += 1; 
        $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'1_1_'.$read}[1]}{b} += 1;    
        print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t$readList{'1_2_'.$read}[0]\t\*\t\*\tBoth";
      }       
      else{print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t$readList{'1_2_'.$read}[0]\t\*\t\*\tError13";}
    } 
  }

  elsif($readList{'1_1_'.$read}) # mate1 aligned to mel but not sec
  {
    if($readList{'2_1_'.$read} && 
       $readList{'2_2_'.$read}) #2
    {
      if($readList{'1_1_'.$read}[0] eq $readList{'2_1_'.$read}[0])
      {
        $geneHash{$readList{'1_1_'.$read}[0]}{$species1} += 1; 
        if($readList{'1_1_'.$read}[1] eq $readList{'2_1_'.$read}[1])
        {
          $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'1_1_'.$read}[1]}{$species1} += 1;
        } 
        else
        {
          $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'1_1_'.$read}[1]}{$species1} += 1; 
          $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'2_1_'.$read}[1]}{$species1} += 1;
        } 
        print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t\*\t$readList{'2_1_'.$read}[0]\t$readList{'2_2_'.$read}[0]\t$species1";
      }
      else{print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t\*\t$readList{'2_1_'.$read}[0]\t$readList{'2_2_'.$read}[0]\tError2";}
    } 
    elsif($readList{'2_1_'.$read}) #3
    {
      if($readList{'1_1_'.$read}[0] eq $readList{'2_1_'.$read}[0])
      {
        $geneHash{$readList{'1_1_'.$read}[0]}{$species1} += 1; 
        if($readList{'1_1_'.$read}[1] eq $readList{'2_1_'.$read}[1])
        {
          $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'1_1_'.$read}[1]}{$species1} += 1;
        } 
        else
        {
          $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'1_1_'.$read}[1]}{$species1} += 1; 
          $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'2_1_'.$read}[1]}{$species1} += 1;
        } 
        print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t\*\t$readList{'2_1_'.$read}[0]\t\*\t$species1";
      }
      else{print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t\*\t$readList{'2_1_'.$read}[0]\t\*\tError3";}
    }
    elsif($readList{'2_2_'.$read}) #14
    {
      print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t\*\t\*\t$readList{'2_2_'.$read}[0]\tTransplicing";
    }
    else #4
    {
      $geneHash{$readList{'1_1_'.$read}[0]}{$species1} += 1; 
      $exonHash{$readList{'1_1_'.$read}[0]}{$readList{'1_1_'.$read}[1]}{$species1} += 1; 
      print OUT1 "\n$read\t$readList{'1_1_'.$read}[0]\t\*\t\*\t\*\t$species1";
    }       
  }

  elsif($readList{'1_2_'.$read}) # mate1 aligned to sec but not mel
  {
    if($readList{'2_1_'.$read} && 
       $readList{'2_2_'.$read}) #7
    {
      if($readList{'1_2_'.$read}[0] eq $readList{'2_2_'.$read}[0])
      {
        $geneHash{$readList{'1_2_'.$read}[0]}{$species2} += 1; 
        if($readList{'1_2_'.$read}[1] eq $readList{'2_2_'.$read}[1])
        {
          $exonHash{$readList{'1_2_'.$read}[0]}{$readList{'1_2_'.$read}[1]}{$species2} += 1;
        } 
        else
        {
          $exonHash{$readList{'1_2_'.$read}[0]}{$readList{'1_2_'.$read}[1]}{$species2} += 1; 
          $exonHash{$readList{'1_2_'.$read}[0]}{$readList{'2_2_'.$read}[1]}{$species2} += 1;
        }  
        print OUT1 "\n$read\t\*\t$readList{'1_2_'.$read}[0]\t$readList{'2_1_'.$read}[0]\t$readList{'2_2_'.$read}[0]\t$species2";
      }
      else{print OUT1 "\n$read\t\*\t$readList{'1_2_'.$read}[0]\t$readList{'2_1_'.$read}[0]\t$readList{'2_2_'.$read}[0]\tError7";}
    } 
    elsif($readList{'2_2_'.$read}) #8
    {
      if($readList{'1_2_'.$read}[0] eq $readList{'2_2_'.$read}[0])
      {
        $geneHash{$readList{'1_2_'.$read}[0]}{$species2} += 1; 
        if($readList{'1_2_'.$read}[1] eq $readList{'2_2_'.$read}[1])
        {
          $exonHash{$readList{'1_2_'.$read}[0]}{$readList{'1_2_'.$read}[1]}{$species2} += 1;
        } 
        else
        {
          $exonHash{$readList{'1_2_'.$read}[0]}{$readList{'1_2_'.$read}[1]}{$species2} += 1; 
          $exonHash{$readList{'1_2_'.$read}[0]}{$readList{'2_2_'.$read}[1]}{$species2} += 1;
        }  
        print OUT1 "\n$read\t\*\t$readList{'1_2_'.$read}[0]\t\*\t$readList{'2_2_'.$read}[0]\t$species2";
      }
      else{print OUT1 "\n$read\t\*\t$readList{'1_2_'.$read}[0]\t\*\t$readList{'2_2_'.$read}[0]\tError8";}
    }
    elsif($readList{$read}{2}{$species1}) #15
    {
      print OUT1 "\n$read\t\*\t$readList{'1_2_'.$read}[0]\t$readList{'2_1_'.$read}[0]\t\*\tTransplicing";
    } 
    else #9
    {
      $geneHash{$readList{'1_2_'.$read}[0]}{$species2} += 1; 
      $exonHash{$readList{'1_2_'.$read}[0]}{$readList{'1_2_'.$read}[1]}{$species2} += 1; 
      print OUT1 "\n$read\t\*\t$readList{'1_2_'.$read}[0]\t\*\t\*\t$species2";        
    }
  }

  else # mate1 aligned to nothing
  {
    if($readList{'2_1_'.$read} && 
       $readList{'2_2_'.$read}) #12
    {
      if($readList{'2_1_'.$read}[0] eq $readList{'2_2_'.$read}[0] && 
         $readList{'2_1_'.$read}[1] eq $readList{'2_2_'.$read}[1])
      {
        $geneHash{$readList{'2_1_'.$read}[0]}{b} += 1; 
        $exonHash{$readList{'2_1_'.$read}[0]}{$readList{'2_1_'.$read}[1]}{b} += 1;   
        print OUT1 "\n$read\t\*\t\*\t$readList{'2_1_'.$read}[0]\t$readList{'2_2_'.$read}[0]\tBoth";
      }
      else{print OUT1 "\n$read\t\*\t\*\t$readList{'2_1_'.$read}[0]\t$readList{'2_2_'.$read}[0]\tError12";}
    }
    elsif($readList{'2_1_'.$read}) #5
    {
      $geneHash{$readList{'2_1_'.$read}[0]}{$species1} += 1; 
      $exonHash{$readList{'2_1_'.$read}[0]}{$readList{'2_1_'.$read}[1]}{$species1} += 1;   
      print OUT1 "\n$read\t\*\t\*\t$readList{'2_1_'.$read}[0]\t\*\t$species1";              
    } 
    elsif($readList{'2_2_'.$read}) #10
    {
      $geneHash{$readList{'2_2_'.$read}[0]}{$species2} += 1; 
      $exonHash{$readList{'2_2_'.$read}[0]}{$readList{'2_2_'.$read}[1]}{$species2} += 1;
      print OUT1 "\n$read\t\*\t\*\t\*\t$readList{'2_2_'.$read}[0]\t$species2";
    }
    else #16
    {
      print OUT1 "\n$read\t\*\t\*\t\*\t\*\tNA";
    }       
  }
}
close OUT1;

open(OUT2,">$out_dir\/$suffix\.mosaik.exons.txt") or die "Error writing to $out_dir\/$suffix\.mosaik.exons.txt\n";
print OUT2 "gene_exon\t$species1\t$species2\tBoth\tAdj$species1\tAdj$species2";

my$madj = 0;
my$sadj = 0;

foreach my$k1 (keys %exonHash)
{
  foreach my$k2 (keys %{$exonHash{$k1}})
  {
    if($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2} == 0)
    {
      print OUT2 "\n$k1\_$k2\t$exonHash{$k1}{$k2}{$species1}\t$exonHash{$k1}{$k2}{$species2}\t$exonHash{$k1}{$k2}{b}\t0\t0";
    }
    else
    {
      $madj = $exonHash{$k1}{$k2}{$species1} + ($exonHash{$k1}{$k2}{$species1} / ($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2})) * $exonHash{$k1}{$k2}{b};  
      $sadj = $exonHash{$k1}{$k2}{$species2} + ($exonHash{$k1}{$k2}{$species2} / ($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2})) * $exonHash{$k1}{$k2}{b};
      print OUT2 "\n$k1\_$k2\t$exonHash{$k1}{$k2}{$species1}\t$exonHash{$k1}{$k2}{$species2}\t$exonHash{$k1}{$k2}{b}\t$madj\t$sadj";
    }
  }
}
close OUT2;

open(OUT3,">$out_dir\/$suffix\.mosaik.genes.txt") or die "Error writing to $out_dir\/$suffix\.mosaik.genes.txt\n";
print OUT3 "gene\t$species1\t$species2\tBoth\tAdj$species1\tAdj$species2";

$madj = 0;
$sadj = 0;

foreach(keys %geneHash)
{
  if($geneHash{$_}{$species1} + $geneHash{$_}{$species2} == 0)
  {
    print OUT3 "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}\t$geneHash{$_}{b}\t0\t0";
  }
  else
  {
    $madj = $geneHash{$_}{$species1} + ($geneHash{$_}{$species1} / ($geneHash{$_}{$species1} + $geneHash{$_}{$species2})) * $geneHash{$_}{b};
    $sadj = $geneHash{$_}{$species2} + ($geneHash{$_}{$species2} / ($geneHash{$_}{$species1} + $geneHash{$_}{$species2})) * $geneHash{$_}{b};
    print OUT3 "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}\t$geneHash{$_}{b}\t$madj\t$sadj";
  }
}
close OUT3;
}
