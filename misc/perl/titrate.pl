#!/usr/bin/perl

########################################################################
# 
# 01/20/2010
#
# titrate.pl
#
# Purpose: sequentially trim reads to improve alignment percentage 
# 
# Input: an initial list of sequence reads 
#
# Output: several alignment files until end of titration is reached, then combine them
#
# Usage: titrate.pl 
#
########################################################################

use strict;
use warnings;
use filetest;

my$fasta1 = $ARGV[0];
my$fasta2 = $ARGV[1];
my$file = $ARGV[2];
my$trim5 = $ARGV[3];
my$trim3 = $ARGV[4];
my$length = $ARGV[5];
my$tool = $ARGV[6];
my$threads = $ARGV[7];
my$st = $ARGV[8];
my$hs = $ARGV[9];

my$ref1;
my$ebwt1;
if($fasta1 =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$ebwt1 = $1;}
if($fasta1 =~ /^(\S*\/){0,}(\S+)\.(fastq|fasta|fa|fq)$/){$ref1 = $2;} 

my$ref2;
my$ebwt2;
if($fasta2 =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$ebwt2 = $1;} 
if($fasta2 =~ /^(\S*\/){0,}(\S+)\.(fastq|fasta|fa|fq)$/){$ref2 = $2;}

my$name;
if($file =~ /^(\S+)\.(fastq|fq|fasta|fa)$/){$name = $1;}

my$cut5; my$cut3; my$current5; my$current3;

# initial alignment to each reference

if($tool == "bowtie")
{

  # trim from the 5' end
  if($trim5 > 0 && $trim3 == 0)
  {
    if($file =~ /^\S+\.(fastq|fq)$/)
    {
      system "./bowtie -q -k 1 -v 0 -p $threads --best --sam $ebwt1 $file --un $name\.$ref1\.bowtie\_trim5_0\.un > $name\.$ref1\.bowtie\_trim5_0\.sam";
      system "./bowtie -q -k 1 -v 0 -p $threads --best --sam $ebwt2 $file --un $name\.$ref2\.bowtie\_trim5_0\.un > $name\.$ref2\.bowtie\_trim5_0\.sam";
    }
    elsif($file =~ /^\S+\.(fasta|fa)$/)
    {
      system "./bowtie -f -k 1 -v 0 -p $threads --best --sam $ebwt1 $file --un $name\.$ref1\.bowtie\_trim5_0\.un > $name\.$ref1\.bowtie\_trim5_0\.sam";
      system "./bowtie -f -k 1 -v 0 -p $threads --best --sam $ebwt2 $file --un $name\.$ref2\.bowtie\_trim5_0\.un > $name\.$ref2\.bowtie\_trim5_0\.sam";
    }

    # titrated alignment
      
    # reference 1
    $cut5 = $trim5;
    $current5 = 0;
    while($length - $cut5 >= 15 && -e "$name\.$ref1\.bowtie\_trim5_$current5\.un") 
    {
      # here is the iterative trimming part
      # 4^14 (as does 4^15) insures uniqueness in Drosophila, maybe needs to change for other sized genomes
      # this is why the smallest allowable read should be >= 15 bases in length
      if($file =~ /^\S+\.(fastq|fq)$/)
      {
        system "./bowtie -q -k 1 -v 0 -p $threads --best --sam --trim5 $cut5 $ebwt1 $name\.$ref1\.bowtie\_trim5_$current5\.un --un $name\.$ref1\.bowtie\_trim5_$cut5\.un > $name\.$ref1\.bowtie\_trim5_$cut5\.sam";
      }
      elsif($file =~ /^\S+\.(fasta|fa)$/)
      {
        system "./bowtie -f -k 1 -v 0 -p $threads --best --sam --trim5 $cut5 $ebwt1 $name\.$ref1\.bowtie\_trim5_$current5\.un --un $name\.$ref1\.bowtie\_trim5_$cut5\.un > $name\.$ref1\.bowtie\_trim5_$cut5\.sam";
      }
      $cut5 += $trim5; 
      $current5 += $trim5;
    }

    # reference 2
    $cut5 = $trim5;
    $current5 = 0;
    while($length - $cut5 >= 15 && -e "$name\.$ref2\.bowtie\_trim5_$current5\.un") 
    {
      # here is the iterative trimming part
      # 4^14 (as does 4^15) insures uniqueness in Drosophila, maybe needs to change for other sized genomes
      # this is why the smallest allowable read should be >= 15 bases in length
      if($file =~ /^\S+\.(fastq|fq)$/)
      {
        system "./bowtie -q -k 1 -v 0 -p $threads --best --sam --trim5 $cut5 $ebwt2 $name\.$ref2\.bowtie\_trim5_$current5\.un --un $name\.$ref2\.bowtie\_trim5_$cut5\.un > $name\.$ref2\.bowtie\_trim5_$cut5\.sam";
      }
      elsif($file =~ /^\S+\.(fasta|fa)$/)
      {
        system "./bowtie -f -k 1 -v 0 -p $threads --best --sam --trim5 $cut5 $ebwt2 $name\.$ref2\.bowtie\_trim5_$current5\.un --un $name\.$ref2\.bowtie\_trim5_$cut5\.un > $name\.$ref2\.bowtie\_trim5_$cut5\.sam";
      }
      $cut5 += $trim5; 
      $current5 += $trim5;
    }
  }

  # trim from the 3' end
  elsif($trim5 == 0 && $trim3 > 0)
  { 
    if($file =~ /^\S+\.(fastq|fq)$/)
    {
      system "./bowtie -q -k 1 -v 0 -p $threads --best --sam $ebwt1 $file --un $name\.$ref1\.bowtie\_trim3_0\.un > $name\.$ref1\.bowtie\_trim3_0\.sam";
      system "./bowtie -q -k 1 -v 0 -p $threads --best --sam $ebwt2 $file --un $name\.$ref2\.bowtie\_trim3_0\.un > $name\.$ref2\.bowtie\_trim3_0\.sam";
    }
    elsif($file =~ /^\S+\.(fasta|fa)$/)
    {
      system "./bowtie -f -k 1 -v 0 -p $threads --best --sam $ebwt1 $file --un $name\.$ref1\.bowtie\_trim3_0\.un > $name\.$ref1\.bowtie\_trim3_0\.sam";
      system "./bowtie -f -k 1 -v 0 -p $threads --best --sam $ebwt2 $file --un $name\.$ref2\.bowtie\_trim3_0\.un > $name\.$ref2\.bowtie\_trim3_0\.sam";
    }

    # titrated alignment
      
    # reference 1
    $cut3 = $trim3;
    $current3 = 0;
    while($length - $cut3 >= 15 && -e "$name\.$ref1\.bowtie\_trim3_$current3\.un") 
    {
      # here is the iterative trimming part
      # 4^14 (as does 4^15) insures uniqueness in Drosophila, maybe needs to change for other sized genomes
      # this is why the smallest allowable read should be >= 15 bases in length
      if($file =~ /^\S+\.(fastq|fq)$/)
      {
        system "./bowtie -q -k 1 -v 0 -p $threads --best --sam --trim3 $cut3 $ebwt1 $name\.$ref1\.bowtie\_trim3_$current3\.un --un $name\.$ref1\.bowtie\_trim3_$cut3\.un > $name\.$ref1\.bowtie\_trim3_$cut3\.sam";
      }
      elsif($file =~ /^\S+\.(fasta|fa)$/)
      {
        system "./bowtie -f -k 1 -v 0 -p $threads --best --sam --trim3 $cut3 $ebwt1 $name\.$ref1\.bowtie\_trim3_$current3\.un --un $name\.$ref1\.bowtie\_trim3_$cut3\.un > $name\.$ref1\.bowtie\_trim3_$cut3\.sam";
      }
      $cut3 += $trim3; 
      $current3 += $trim3;
    }

    # reference 2
    $cut3 = $trim3;
    $current3 = 0;
    while($length - $cut3 >= 15 && -e "$name\.$ref1\.bowtie\_trim3_$current3\.un") 
    {
      # here is the iterative trimming part
      # 4^14 (as does 4^15) insures uniqueness in Drosophila, maybe needs to change for other sized genomes
      # this is why the smallest allowable read should be >= 15 bases in length
      if($file =~ /^\S+\.(fastq|fq)$/)
      {
        system "./bowtie -q -k 1 -v 0 -p $threads --best --sam --trim3 $cut3 $ebwt2 $name\.$ref2\.bowtie\_trim3_$current3\.un --un $name\.$ref2\.bowtie\_trim3_$cut3\.un > $name\.$ref2\.bowtie\_trim3_$cut3\.sam";
      }
      elsif($file =~ /^\S+\.(fasta|fa)$/)
      {
        system "./bowtie -f -k 1 -v 0 -p $threads --best --sam --trim3 $cut3 $ebwt2 $name\.$ref2\.bowtie\_trim3_$current3\.un --un $name\.$ref2\.bowtie\_trim3_$cut3\.un > $name\.$ref2\.bowtie\_trim3_$cut3\.sam";
      }
      $cut3 += $trim3; 
      $current3 += $trim3;
    }
  }

  # trim from the 5' and 3' ends
  elsif($trim5 > 0 && $trim3 > 0)
  { 
if($file =~ /^\S+\.(fastq|fq)$/)
    {
      system "./bowtie -q -k 1 -v 0 -p $threads --best --sam $ebwt1 $file --un $name\.$ref1\.bowtie\_trim5_0_trim3_0\.un > $name\.$ref1\.bowtie\_trim5_0_trim3_0\.sam";
      system "./bowtie -q -k 1 -v 0 -p $threads --best --sam $ebwt2 $file --un $name\.$ref2\.bowtie\_trim5_0_trim3_0\.un > $name\.$ref2\.bowtie\_trim5_0_trim3_0\.sam";
    }
    elsif($file =~ /^\S+\.(fasta|fa)$/)
    {
      system "./bowtie -f -k 1 -v 0 -p $threads --best --sam $ebwt1 $file --un $name\.$ref1\.bowtie\_trim5_0_trim3_0\.un > $name\.$ref1\.bowtie\_trim5_0_trim3_0\.sam";
      system "./bowtie -f -k 1 -v 0 -p $threads --best --sam $ebwt2 $file --un $name\.$ref2\.bowtie\_trim5_0_trim3_0\.un > $name\.$ref2\.bowtie\_trim5_0_trim3_0\.sam";
    }

    # titrated alignment
      
    # reference 1
    $cut5 = $trim5;
    $current5 = 0;
    $cut3 = $trim3;
    $current3 = 0;
    while($length - $cut5 - $cut3 >= 15 && -e "$name\.$ref1\.bowtie\_trim5_$current5\_trim3_$current3\.un") 
    {
      # here is the iterative trimming part
      # 4^14 (as does 4^15) insures uniqueness in Drosophila, maybe needs to change for other sized genomes
      # this is why the smallest allowable read should be >= 15 bases in length
      if($file =~ /^\S+\.(fastq|fq)$/)
      {
        system "./bowtie -q -k 1 -v 0 -p $threads --best --sam --trim5 $cut5 --trim3 $cut3 $ebwt1 $name\.$ref1\.bowtie\_trim5_$current5\_trim3_$current3\.un --un $name\.$ref1\.bowtie\_trim5_$cut5\_trim3_$cut3\.un > $name\.$ref1\.bowtie\_trim5_$cut5\_trim3_$cut3\.sam";
      }
      elsif($file =~ /^\S+\.(fasta|fa)$/)
      {
        system "./bowtie -f -k 1 -v 0 -p $threads --best --sam --trim5 $cut5 --trim3 $cut3 $ebwt1 $name\.$ref1\.bowtie\_trim5_$current5\_trim3_$current3\.un --un $name\.$ref1\.bowtie\_trim5_$cut5\_trim3_$cut3\.un > $name\.$ref1\.bowtie\_trim5_$cut5\_trim3_$cut3\.sam";
      }
      $cut5 += $trim5; 
      $current5 += $trim5;
      $cut3 += $trim3;
      $current3 += $trim3;
    }

    # reference 2
    $cut5 = $trim5;
    $current5 = 0;
    $cut3 = $trim3;
    $current3 = 0;
    while($length - $cut5 - $cut3 >= 15 && -e "$name\.$ref2\.bowtie\_trim5_$current5\_trim3_$current3\.un") 
    {
      # here is the iterative trimming part
      # 4^14 (as does 4^15) insures uniqueness in Drosophila, maybe needs to change for other sized genomes
      # this is why the smallest allowable read should be >= 15 bases in length
      if($file =~ /^\S+\.(fastq|fq)$/)
      {
        system "./bowtie -q -k 1 -v 0 -p $threads --best --sam --trim5 $cut5 --trim3 $cut3 $ebwt2 $name\.$ref2\.bowtie\_trim5_$current5\_trim3_$current3\.un --un $name\.$ref2\.bowtie\_trim5_$cut5\_trim3_$cut3\.un > $name\.$ref2\.bowtie\_trim5_$cut5\_trim3_$cut3\.sam";
      }
      elsif($file =~ /^\S+\.(fasta|fa)$/)
      {
        system "./bowtie -f -k 1 -v 0 -p $threads --best --sam --trim5 $cut5 --trim3 $cut3 $ebwt2 $name\.$ref2\.bowtie\_trim5_$current5\_trim3_$current3\.un --un $name\.$ref2\.bowtie\_trim5_$cut5\_trim3_$cut3\.un > $name\.$ref2\.bowtie\_trim5_$cut5\_trim3_$cut3\.sam";
      }
      $cut5 += $trim5; 
      $current5 += $trim5;
      $cut3 += $trim3;
      $current3 += $trim3;
    }
  }
}
elsif($tool == "mosaik")
{
  # trim from the 5' end
  if($trim5 > 0 && $trim3 == 0)
  {
    print "\nBuilding $file...\n";
    if($file =~ /^\S+\.(fastq|fq)$/)
    {
      system "MosaikBuild -q $file -out $name\.dat -st $st";
      system "MosaikAligner -in $name\.dat -out $name\.$ref1\.mosaik\_trim5_0\.dat -rur $name\.$ref1\.trim5_0\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\_trim5_0\.output";
      system "MosaikAligner -in $name\.dat -out $name\.$ref2\.mosaik\_trim5_0\.dat -rur $name\.$ref2\.trim5_0\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\_trim5_0\.output";
      system "MosaikText -in $name\.$ref1\.mosaik\_trim5_0\.dat -bed $name\.$ref1\.mosaik\_trim5_0\.bed";
      system "MosaikText -in $name\.$ref2\.mosaik\_trim5_0\.dat -bed $name\.$ref2\.mosaik\_trim5_0\.bed";
    }
    elsif($file =~ /^\S+\.(fasta|fa)$/)
    {
      system "MosaikBuild -fr $file -out $name\.dat -st illumina -assignQual 60";
      system "MosaikAligner -in $name\.dat -out $name\.$ref1\.mosaik\_trim5_0\.dat -rur $name\.$ref1\.trim5_0\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\_trim5_0\.output";
      system "MosaikAligner -in $name\.dat -out $name\.$ref2\.mosaik\_trim5_0\.dat -rur $name\.$ref2\.trim5_0\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\_trim5_0\.output";
      system "MosaikText -in $name\.$ref1\.mosaik\_trim5_0\.dat -bed $name\.$ref1\.mosaik\_trim5_0\.bed";
      system "MosaikText -in $name\.$ref2\.mosaik\_trim5_0\.dat -bed $name\.$ref2\.mosaik\_trim5_0\.bed";
    }

    # titrated alignment
      
    # reference 1
    $cut5 = $trim5;
    $current5 = 0;
    while($length - $cut5 >= 15 && -e "$name\.$ref1\.trim5_$current5\.un") 
    {
      # here is the iterative trimming part
      # 4^14 (as does 4^15) insures uniqueness in Drosophila, maybe needs to change for other sized genomes
      # this is why the smallest allowable read should be >= 15 bases in length

      system "perl read_trimmer.pl $name\.$ref1\.trim5_$current5\.un $cut5 0 $name\.$ref1\.trim5_$cut5\.fq";

      system "MosaikBuild -q $name\.$ref1\.trim5_$cut5\.fq -out $name\.$ref1\.trim5_$cut5\.dat -st $st";
      system "MosaikAligner -in $name\.$ref1\.trim5_$cut5\.dat -out $name\.$ref1\.mosaik\_trim5_$cut5\.dat -rur $name\.$ref1\.trim5_$cut5\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\_trim5_$cut5\.output";
      system "MosaikText -in $name\.$ref1\.mosaik\_trim5_$cut5\.dat -bed $name\.$ref1\.mosaik\_trim5_$cut5\.bed";

      $cut5 += $trim5; 
      $current5 += $trim5;
    }

    # reference 2
    $cut5 = $trim5;
    $current5 = 0;
    while($length - $cut5 >= 15 && -e "$name\.$ref2\.bowtie\_trim5_$current5\.un") 
    {
      # here is the iterative trimming part
      # 4^14 (as does 4^15) insures uniqueness in Drosophila, maybe needs to change for other sized genomes
      # this is why the smallest allowable read should be >= 15 bases in length

      system "perl read_trimmer.pl $name\.$ref2\.trim5_$current5\.un $cut5 0 $name\.$ref2\.trim5_$cut5\.fq";

      system "MosaikBuild -q $name\.$ref2\.trim5_$cut5\.fq -out $name\.$ref2\.trim5_$cut5\.dat -st $st";
      system "MosaikAligner -in $name\.$ref2\.trim5_$cut5\.dat -out $name\.$ref2\.mosaik\_trim5_$cut5\.dat -rur $name\.$ref2\.trim5_$cut5\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\_trim5_$cut5\.output";
      system "MosaikText -in $name\.$ref2\.mosaik\_trim5_$cut5\.dat -bed $name\.$ref2\.mosaik\_trim5_$cut5\.bed";

      $cut5 += $trim5; 
      $current5 += $trim5;
    }
  }
  elsif($trim5 == 0 && $trim3 > 0)
  {
    print "\nBuilding $file...\n";
    if($file =~ /^\S+\.(fastq|fq)$/)
    {
      system "MosaikBuild -q $file -out $name\.dat -st $st";
      system "MosaikAligner -in $name\.dat -out $name\.$ref1\.mosaik\_trim3_0\.dat -rur $name\.$ref1\.trim3_0\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\_trim3_0\.output";
      system "MosaikAligner -in $name\.dat -out $name\.$ref2\.mosaik\_trim3_0\.dat -rur $name\.$ref2\.trim3_0\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\_trim3_0\.output";
      system "MosaikText -in $name\.$ref1\.mosaik\_trim3_0\.dat -bed $name\.$ref1\.mosaik\_trim3_0\.bed";
      system "MosaikText -in $name\.$ref2\.mosaik\_trim3_0\.dat -bed $name\.$ref2\.mosaik\_trim3_0\.bed";
    }
    elsif($file =~ /^\S+\.(fasta|fa)$/)
    {
      system "MosaikBuild -fr $file -out $name\.dat -st illumina -assignQual 60";
      system "MosaikAligner -in $name\.dat -out $name\.$ref1\.mosaik\_trim3_0\.dat -rur $name\.$ref1\.trim3_0\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\_trim3_0\.output";
      system "MosaikAligner -in $name\.dat -out $name\.$ref2\.mosaik\_trim3_0\.dat -rur $name\.$ref2\.trim3_0\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\_trim3_0\.output";
      system "MosaikText -in $name\.$ref1\.mosaik\_trim3_0\.dat -bed $name\.$ref1\.mosaik\_trim3_0\.bed";
      system "MosaikText -in $name\.$ref2\.mosaik\_trim3_0\.dat -bed $name\.$ref2\.mosaik\_trim3_0\.bed";
    }

    # titrated alignment
      
    # reference 1
    $cut3 = $trim3;
    $current3 = 0;
    while($length - $cut3 >= 15 && -e "$name\.$ref1\.trim3_$current3\.un") 
    {
      # here is the iterative trimming part
      # 4^14 (as does 4^15) insures uniqueness in Drosophila, maybe needs to change for other sized genomes
      # this is why the smallest allowable read should be >= 15 bases in length

      system "perl read_trimmer.pl $name\.$ref1\.trim3_$current3\.un 0 $cut3 $name\.$ref1\.trim3_$cut3\.fq";

      system "MosaikBuild -q $name\.$ref1\.trim3_$cut3\.fq -out $name\.$ref1\.trim3_$cut3\.dat -st $st";
      system "MosaikAligner -in $name\.$ref1\.trim3_$cut3\.dat -out $name\.$ref1\.mosaik\_trim3_$cut3\.dat -rur $name\.$ref1\.trim3_$cut3\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\_trim3_$cut3\.output";
      system "MosaikText -in $name\.$ref1\.mosaik\_trim3_$cut3\.dat -bed $name\.$ref1\.mosaik\_trim3_$cut3\.bed";

      $cut3 += $trim3; 
      $current3 += $trim3;
    }

    # reference 2
    $cut3 = $trim3;
    $current3 = 0;
    while($length - $cut3 >= 15 && -e "$name\.$ref2\.bowtie\_trim3_$current3\.un") 
    {
      # here is the iterative trimming part
      # 4^14 (as does 4^15) insures uniqueness in Drosophila, maybe needs to change for other sized genomes
      # this is why the smallest allowable read should be >= 15 bases in length

      system "perl read_trimmer.pl $name\.$ref2\.trim3_$current3\.un 0 $cut3 $name\.$ref2\.trim3_$cut3\.fq";

      system "MosaikBuild -q $name\.$ref2\.trim3_$cut3\.fq -out $name\.$ref2\.trim3_$cut3\.dat -st $st";
      system "MosaikAligner -in $name\.$ref2\.trim3_$cut3\.dat -out $name\.$ref2\.mosaik\_trim3_$cut3\.dat -rur $name\.$ref2\.trim3_$cut3\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\_trim3_$cut3\.output";
      system "MosaikText -in $name\.$ref2\.mosaik\_trim3_$cut3\.dat -bed $name\.$ref2\.mosaik\_trim3_$cut3\.bed";

      $cut3 += $trim3; 
      $current3 += $trim3;
    }
  }
  elsif($trim5 > 0 && $trim3 > 0)
  {
    print "\nBuilding $file...\n";
    if($file =~ /^\S+\.(fastq|fq)$/) 
    {
      system "MosaikBuild -q $file -out $name\.dat -st $st";
      system "MosaikAligner -in $name\.dat -out $name\.$ref1\.mosaik\_trim5_0_trim3_0\.dat -rur $name\.$ref1\.trim5_0_trim3_0\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\_trim5_0_trim3_0\.output";
      system "MosaikAligner -in $name\.dat -out $name\.$ref2\.mosaik\_trim5_0_trim3_0\.dat -rur $name\.$ref2\.trim5_0_trim3_0\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\_trim5_0_trim3_0\.output";
      system "MosaikText -in $name\.$ref1\.mosaik\_trim5_0_trim3_0\.dat -bed $name\.$ref1\.mosaik\_trim5_0_trim3_0\.bed";
      system "MosaikText -in $name\.$ref2\.mosaik\_trim5_0_trim3_0\.dat -bed $name\.$ref2\.mosaik\_trim5_0_trim3_0\.bed";
    }
    elsif($file =~ /^\S+\.(fasta|fa)$/)
    {
      system "MosaikBuild -fr $file -out $name\.dat -st illumina -assignQual 60";
      system "MosaikAligner -in $name\.dat -out $name\.$ref1\.mosaik\_trim5_0_trim3_0\.dat -rur $name\.$ref1\.trim5_0_trim3_0\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\_trim5_0_trim3_0\.output";
      system "MosaikAligner -in $name\.dat -out $name\.$ref2\.mosaik\_trim5_0_trim3_0\.dat -rur $name\.$ref2\.trim5_0_trim3_0\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\_trim5_0_trim3_0\.output";
      system "MosaikText -in $name\.$ref1\.mosaik\_trim5_0_trim3_0\.dat -bed $name\.$ref1\.mosaik\_trim5_0_trim3_0\.bed";
      system "MosaikText -in $name\.$ref2\.mosaik\_trim5_0_trim3_0\.dat -bed $name\.$ref2\.mosaik\_trim5_0_trim3_0\.bed";
    }

    # titrated alignment
      
    # reference 1
    $cut5 = $trim5;
    $current5 = 0;
    $cut3 = $trim3;
    $current3 = 0;
    while($length - $cut5 - $cut3 >= 15 && -e "$name\.$ref1\.trim5_$current5\_trim3_$current3\.un") 
    {
      # here is the iterative trimming part
      # 4^14 (as does 4^15) insures uniqueness in Drosophila, maybe needs to change for other sized genomes
      # this is why the smallest allowable read should be >= 15 bases in length

      system "perl read_trimmer.pl $name\.$ref1\.trim5_$current5\_trim3_$current3\.un $cut5 0 $name\.$ref1\.trim5_$cut5\_trim3_$cut3\.fq";

      system "MosaikBuild -q $name\.$ref1\.trim5_$cut5\_trim3_$cut3\.fq -out $name\.$ref1\.trim5_$cut5\_trim3_$cut3\.dat -st $st";
      system "MosaikAligner -in $name\.$ref1\.trim5_$cut5\_trim3_$cut3\.dat -out $name\.$ref1\.mosaik\_trim5_$cut5\_trim3_$cut3\.dat -rur $name\.$ref1\.trim5_$cut5\_trim3_$cut3\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\_trim5_$cut5\_trim3_$cut3\.output";
      system "MosaikText -in $name\.$ref1\.mosaik\_trim5_$cut5\_trim3_$cut3\.dat -bed $name\.$ref1\.mosaik\_trim5_$cut5\_trim3_$cut3\.bed";

      $cut5 += $trim5; 
      $current5 += $trim5;
      $cut3 += $trim3;
      $current3 += $trim3;
    }

    # reference 2
    $cut5 = $trim5;
    $current5 = 0;
    $cut3 = $trim3;
    $current3 = 0;
    while($length - $cut5 - $cut3 >= 15 && -e "$name\.$ref2\.bowtie\_trim5_$current5\_trim3_$current3\.un") 
    {
      # here is the iterative trimming part
      # 4^14 (as does 4^15) insures uniqueness in Drosophila, maybe needs to change for other sized genomes
      # this is why the smallest allowable read should be >= 15 bases in length

      system "perl read_trimmer.pl $name\.$ref2\.trim5_$current5\_trim3_$current3\.un $cut5 $cut3 $name\.$ref2\.trim5_$cut5\_trim3_$cut3\.fq";

      system "MosaikBuild -q $name\.$ref2\.trim5_$cut5\_trim3_$cut3\.fq -out $name\.$ref2\.trim5_$cut5\_trim3_$cut3\.dat -st $st";
      system "MosaikAligner -in $name\.$ref2\.trim5_$cut5\_trim3_$cut3\.dat -out $name\.$ref2\.mosaik\_trim5_$cut5\_trim3_$cut3\.dat -rur $name\.$ref2\.trim5_$cut5\_trim3_$cut3\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\_trim5_$cut5\_trim3_$cut3\.output";
      system "MosaikText -in $name\.$ref2\.mosaik\_trim5_$cut5\_trim3_$cut3\.dat -bed $name\.$ref2\.mosaik\_trim5_$cut5\_trim3_$cut3\.bed";

      $cut5 += $trim5; 
      $current5 += $trim5;
      $cut3 += $trim3;
      $current3 += $trim3;
    }
  }
}

