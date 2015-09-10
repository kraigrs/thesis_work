#!/usr/bin/perl

########################################################################
# 
# 03/31/2011
#
# RegDiv.pl
#
# Purpose: elucidate regulatory divergence 
#          within and between species by quantifying allele-specific 
#          expression from next-generation sequencing data 
# 
# Input: [options]* 
#        -ref1 <ref1> -ref2 <ref2> 
#        {-mix <mix> | -mix1 <mix1> -mix2 <mix2>} 
#        {-hyb1 <hyb1> | -hyb1_1 <hyb1_1> -hyb1_2 <hyb1_2>} 
#        {-hyb2 <hyb2> | -hyb2_1 <hyb2_1> -hyb2_2 <hyb2_2>}
#
# Output: allele-specific expression values for mixed parents and hybrids
#
# Notes: This is a revamping of the original "GUI" pipeline intended to be run in a single command
#
# Fixes: 
#
########################################################################

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

main();
sub main
{

# print error messages, help, and manual
my $man = '';
my $help = '';

# declare defaults
# essential
my$fasta1; my$fasta2;
my$mix; my$mix1; my$mix2;
my$hyb1; my$hyb1_1; my$hyb1_2;
my$hyb2; my$hyb2_1; my$hyb2_2;
my@files2align;

# non-essential
my$form;
my$tool = 'bowtie'; # alignment tool
my$threads = 1;     # num processors to run
my$length = '';     # length of reads
my$st = 'illumina'; # sequencing technology
my$hs = 14;         # hash size (14 guarantees uniqueness in drosophila genome)
my$titrate = '';    # sequentially trim reads
my$short = '';      # trim reads 
my$trim3 = 0;       # num bases to trim from 3' end
my$trim5 = 0;       # num bases to trim from 5' end

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested

if(@ARGV > 0)
{
  GetOptions('help|?'    => \$help, 
             'man'       => \$man,
             'form=s'    => \$form,
             'ref1=s'    => \$fasta1,
             'ref2=s'    => \$fasta2,
             'mix=s'     => \$mix,
             'mix1=s'    => \$mix1,
             'mix2=s'    => \$mix2,
             'hyb1=s'    => \$hyb1,
             'hyb1_1=s'  => \$hyb1_1,
             'hyb1_2=s'  => \$hyb1_2,
             'hyb2=s'    => \$hyb2,
             'hyb2_1=s'  => \$hyb2_1,
             'hyb2_2=s'  => \$hyb2_2,
             'tool=s'    => \$tool,
             'threads=i' => \$threads,
             'titrate'   => \$titrate,
             'short'     => \$short,
             'length=i'  => \$length,
             'trim3=i'   => \$trim3,
             'trim5=i'   => \$trim5,
             'st=s'      => \$st,
             'hs=i'      => \$hs) 
  or pod2usage(-verbose => 1, -input => "RegDiv.pod");
}
else{pod2usage(-verbose => 1, -input => "RegDiv.pod");}

if($help or $man) 
{
  pod2usage(-verbose => 1, -input => "RegDiv.pod") if $help;
  pod2usage(-verbose => 2, -input => "RegDiv.pod") if $man;
}

# specify errors for missing essential values
unless($fasta1 && $fasta2){}

my$ref1;
my$ebwt1;
if($fasta1 =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$ebwt1 = $1;}
if($fasta1 =~ /^(\S*\/){0,}(\S+)\.(fastq|fasta|fa|fq)$/){$ref1 = $2;} 
#print "\nFasta: $fasta1\tRef: $ref1\tEBWT: $ebwt1\n";

my$ref2;
my$ebwt2;
if($fasta2 =~ /^(\S+)\.(fastq|fasta|fa|fq)$/){$ebwt2 = $1;} 
if($fasta2 =~ /^(\S*\/){0,}(\S+)\.(fastq|fasta|fa|fq)$/){$ref2 = $2;} 
#print "\nFasta: $fasta2\tRef: $ref2\tEBWT: $ebwt2\n";

# alignment strategy

my$name;
my$file;

if($tool eq "bowtie")
{
  if($titrate)
  {
    print "Using Bowtie to index $ebwt1 and $ebwt2\n";
    system "bowtie-build $fasta1 $ebwt1";
    system "bowtie-build $fasta2 $ebwt2";

    foreach(@files2align)
    {
      $file = $_;
      system "perl titrate.pl $fasta1 $fasta2 $file $trim5 $trim3 $length $tool $threads null null";
    }
  }
  else
  {
    print "Using Bowtie to index $ebwt1 and $ebwt2\n";
    system "bowtie-build $fasta1 $ebwt1";
    system "bowtie-build $fasta2 $ebwt2";
   
    if($short)
    {
      print "\nAligning $file to $ref1 and $ref2...\n";
      foreach(@files2align)
      {
        $file = $_;
        if($file =~ /^(\S+)\.(fastq|fq|fasta|fa)$/){$name = $1;}
        if($file =~ /^\S+\.(fastq|fq)$/)
        {
          system "./bowtie -q -k 1 -v 0 -p $threads --best --sam --trim5 $trim5 --trim3 $trim3 $ebwt1 $file --un $name\.$ref1\.bowtie\_trim5_$trim5\_trim3_$trim3\.un > $name\.$ref1\.bowtie\_trim5_$trim5\_trim3_$trim3\.sam";
          system "./bowtie -q -k 1 -v 0 -p $threads --best --sam --trim5 $trim5 --trim3 $trim3 $ebwt2 $file --un $name\.$ref2\.bowtie\_trim5_$trim5\_trim3_$trim3\.un > $name\.$ref2\.bowtie\_trim5_$trim5\_trim3_$trim3\.sam";
        }
        elsif($file =~ /^\S+\.(fasta|fa)$/)
        {
          system "./bowtie -f -k 1 -v 0 -p $threads --best --sam --trim5 $trim5 --trim3 $trim3 $ebwt1 $file --un $name\.$ref1\.bowtie\_trim5_$trim5\_trim3_$trim3\.un > $name\.$ref1\.bowtie\_trim5_$trim5\_trim3_$trim3\.sam";
          system "./bowtie -f -k 1 -v 0 -p $threads --best --sam --trim5 $trim5 --trim3 $trim3 $ebwt2 $file --un $name\.$ref2\.bowtie\_trim5_$trim5\_trim3_$trim3\.un > $name\.$ref2\.bowtie\_trim5_$trim5\_trim3_$trim3\.sam";
        }
      }
    }
    else
    {
      print "Using Bowtie to index $ebwt1 and $ebwt2\n";
      system "bowtie-build $fasta1 $ebwt1";
      system "bowtie-build $fasta2 $ebwt2";

      print "\nAligning $file to $ref1 and $ref2...\n";
      foreach(@files2align)
      {
        $file = $_;
        if($file =~ /^(\S+)\.(fastq|fq|fasta|fa)$/){$name = $1;}

        if($file =~ /^\S+\.(fastq|fq)$/)
        {
          system "./bowtie -q -k 1 -v 0 -p $threads --best --sam $ebwt1 $file --un $name\.$ref1\.bowtie\.un > $name\.$ref1\.bowtie\.sam";
          system "./bowtie -q -k 1 -v 0 -p $threads --best --sam $ebwt2 $file --un $name\.$ref2\.bowtie\.un > $name\.$ref2\.bowtie\.sam";
        }
        elsif($file =~ /^\S+\.(fasta|fa)$/)
        {
          system "./bowtie -f -k 1 -v 0 -p $threads --best --sam $ebwt1 $file --un $name\.$ref1\.bowtie\.un > $name\.$ref1\.bowtie\.sam";
          system "./bowtie -f -k 1 -v 0 -p $threads --best --sam $ebwt2 $file --un $name\.$ref2\.bowtie\.un > $name\.$ref2\.bowtie\.sam";
        }
      }
    }
  }
}
elsif($tool eq "mosaik")
{  
  if($titrate)
  {
    print "Using MOSAIK to build reference genomes $ref1 and $ref2\n";
    system "MosaikBuild -fr $fasta1 -oa $ebwt1\.dat";
    system "MosaikBuild -fr $fasta2 -oa $ebwt2\.dat";
    system "MosaikJump -ia $ebwt1\.dat -out $ebwt1\_$hs -hs $hs";
    system "MosaikJump -ia $ebwt2\.dat -out $ebwt2\_$hs -hs $hs";

    foreach(@files2align)
    {
      $file = $_;
      system "perl titrate.pl $fasta1 $fasta2 $file $trim5 $trim3 $length $tool $threads $st $hs";
    }
  }
  else
  {
    if($short)
    {
      print "Using MOSAIK to build reference genomes $ref1 and $ref2\n";
      system "MosaikBuild -fr $fasta1 -oa $ebwt1\.dat";
      system "MosaikBuild -fr $fasta2 -oa $ebwt2\.dat";
      system "MosaikJump -ia $ebwt1\.dat -out $ebwt1\_$hs -hs $hs";
      system "MosaikJump -ia $ebwt2\.dat -out $ebwt2\_$hs -hs $hs";

      foreach(@files2align)
      {
        $file = $_;
        if($file =~ /^(\S+)\.(fastq|fq|fasta|fa)$/){$name = $1;}
        print "\nBuilding and trimming $file...\n";
        if($file =~ /^\S+\.(fastq|fq)$/)
        {
          system "perl read_trimmer.pl $file $trim5 $trim3 $name\.trim5_$trim5\_trim3_$trim3\.fq";
          system "MosaikBuild -q $name\.$ref1\.trim5_$trim5\_trim3_$trim3\.fq -out $name\.$ref1\.trim5_$trim5\_trim3_$trim3\.dat -st $st";

          print "\nAligning $file to $ref1 and $ref2...\n";
          system "MosaikAligner -in $name\.trim5_$trim5\_trim3_$trim3\.dat -out $name\.$ref1\.mosaik\_trim5_$trim5\_trim3_$trim3\.dat -rur $name\.$ref1\.trim5_$trim5\_trim3_$trim3\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\_trim5_$trim5\_trim3_$trim3\.output";
          system "MosaikAligner -in $name\.trim5_$trim5\_trim3_$trim3\.dat -out $name\.$ref2\.mosaik\_trim5_$trim5\_trim3_$trim3\.dat -rur $name\.$ref2\.trim5_$trim5\_trim3_$trim3\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\_trim5_$trim5\_trim3_$trim3\.output";

          print "\nConverting $file alignments into BED format...\n";
          system "MosaikText -in $name\.$ref1\.mosaik\_trim5_$trim5\_trim3_$trim3\.dat -bed $name\.$ref1\.mosaik\_trim5_$trim5\_trim3_$trim3\.bed";
          system "MosaikText -in $name\.$ref2\.mosaik\_trim5_$trim5\_trim3_$trim3\.dat -bed $name\.$ref2\.mosaik\_trim5_$trim5\_trim3_$trim3\.bed";
        }
        elsif($file =~ /^\S+\.(fasta|fa)$/)
        {
          system "perl read_trimmer.pl $file $trim5 $trim3 $name\.trim5_$trim5\_trim3_$trim3\.fa";
          system "MosaikBuild -fr $name\.trim5_$trim5\_trim3_$trim3\.fa -out $name\.trim5_$trim5\_trim3_$trim3\.dat -st illumina -assignQual 60";

          print "\nAligning $file to $ref1 and $ref2...\n";
          system "MosaikAligner -in $name\.trim5_$trim5\_trim3_$trim3\.dat -out $name\.$ref1\.mosaik\_trim5_$trim5\_trim3_$trim3\.dat -rur $name\.$ref1\.trim5_$trim5\_trim3_$trim3\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\_trim5_$trim5\_trim3_$trim3\.output";
          system "MosaikAligner -in $name\.trim5_$trim5\_trim3_$trim3\.dat -out $name\.$ref2\.mosaik\_trim5_$trim5\_trim3_$trim3\.dat -rur $name\.$ref2\.trim5_$trim5\_trim3_$trim3\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\_trim5_$trim5\_trim3_$trim3\.output";

          print "\nConverting $file alignments into BED format...\n";
          system "MosaikText -in $name\.$ref1\.mosaik\_trim5_$trim5\_trim3_$trim3\.dat -bed $name\.$ref1\.mosaik\_trim5_$trim5\_trim3_$trim3\.bed";
          system "MosaikText -in $name\.$ref2\.mosaik\_trim5_$trim5\_trim3_$trim3\.dat -bed $name\.$ref2\.mosaik\_trim5_$trim5\_trim3_$trim3\.bed";
        }
      }
    }
    else
    {
      print "Using MOSAIK to build reference genomes $ref1 and $ref2\n";
      system "MosaikBuild -fr $fasta1 -oa $ebwt1\.dat";
      system "MosaikBuild -fr $fasta2 -oa $ebwt2\.dat";
      system "MosaikJump -ia $ebwt1\.dat -out $ebwt1\_$hs -hs $hs";
      system "MosaikJump -ia $ebwt2\.dat -out $ebwt2\_$hs -hs $hs";

      foreach(@files2align)
      {
        $file = $_;
        if($file =~ /^(\S+)\.(fastq|fq|fasta|fa)$/){$name = $1;}
        if($file =~ /^\S+\.(fastq|fq)$/)
        {
          print "\nBuilding $file...\n";
          system "MosaikBuild -q $file -out $name\.dat -st $st";

          print "\nAligning $file to $ref1 and $ref2...\n";
          system "MosaikAligner -in $name\.dat -out $name\.$ref1\.mosaik\.dat -rur $name\.$ref1\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\.output";
          system "MosaikAligner -in $name\.dat -out $name\.$ref2\.mosaik\.dat -rur $name\.$ref2\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\.output";

          print "\nConverting $file alignments into BED format...\n";
          system "MosaikText -in $name\.$ref1\.mosaik\.dat -bed $name\.$ref1\.mosaik\.bed";
          system "MosaikText -in $name\.$ref2\.mosaik\.dat -bed $name\.$ref2\.mosaik\.bed";
        }
        elsif($file =~ /^\S+\.(fasta|fa)$/)
        {
          print "\nBuilding $file...\n";
          system "MosaikBuild -fr $file -out $name\.dat -st illumina -assignQual 60";

          print "\nAligning $file to $ref1 and $ref2...\n";
          system "MosaikAligner -in $name\.dat -out $name\.$ref1\.mosaik\.dat -rur $name\.$ref1\.un -ia $ebwt1\.dat -hs $hs -mm 0 -m unique -j $ebwt1\_$hs -p $threads > $name\.$ref1\.mosaik\.output";
          system "MosaikAligner -in $name\.dat -out $name\.$ref2\.mosaik\.dat -rur $name\.$ref2\.un -ia $ebwt2\.dat -hs $hs -mm 0 -m unique -j $ebwt2\_$hs -p $threads > $name\.$ref2\.mosaik\.output";

          print "\nConverting $file alignments into BED format...\n";
          system "MosaikText -in $name\.$ref1\.mosaik\.dat -bed $name\.$ref1\.mosaik\.bed";
          system "MosaikText -in $name\.$ref2\.mosaik\.dat -bed $name\.$ref2\.mosaik\.bed";
        }
      }
    }
  }
}

}

