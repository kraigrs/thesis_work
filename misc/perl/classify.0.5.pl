#!/usr/bin/perl

###################################################################################
#
# 12/03/2010
#
# classify.0.5.pl
#
# Purpose: combine 4 .bed files (lifted) to elucidate allele-specific expression 
# 
# Input: two separate alignment files (1 for each species)
#
# Output: a list of genes and their respective match and mismatch counts, along with a 
#         summary of all the reads 
#
# Fixes: this script represents an updated version fixing a memory leak problem. Also, instead of repeatedly
#        building gene and exon hashes, those are made once in the beginning. The memory problem was fixed using 
#        the DB_File module 
#
# E.g. perl classify.0.2.pl <const> <species1> <species2> <bedmate1species1> <bedmate2species1> <bedmate1species2> <bedmate2species2> <readout> <exonout> <geneout>
#     
#      <const>            ==> list of constitutive exons
#      <species1>         ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <species2>         ==> indicate name of 2nd species aligned to
#      <bedmate1species1> ==> file with mate 1 alignments to 1st species
#      <bedmate2species1> ==> file with mate 2 alignments to 1st species
#      <bedmate1species2> ==> file with mate 1 alignments to 2nd species
#      <bedmate2species2> ==> file with mate 2 alignments to 2nd species
#      <readout>          ==> file to send reads for output
#      <exonout>          ==> file to send exon counts for output
#      <geneout>          ==> file to send gene counts for output
#
# perl classify.0.5.pl ../McManus/constitutive_regions_overlap_filtered.bed zhr z30 ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.zhr.mosaik.bed.lifted.converted ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.zhr.mosaik.bed.lifted.converted ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate1.z30.mosaik.bed.lifted.converted ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mate2.z30.mosaik.bed.lifted.converted ../mel_mel_data/mRNA-Seq/z30Xzhr z30Xzhr_v2 &
#
###################################################################################

use strict;
use warnings;
use DB_File;
use Fcntl;

main();
sub main
{

our($line,$value,$gene,$exon,$read,$name);
my@elements;
our(%geneHash,%exonHash,%readNames,%readList,%readList11,%readList12,%readList21,%readList22);

my$const = $ARGV[0];
my$species1 = $ARGV[1]; # Dmel
my$species2 = $ARGV[2]; # Dsec
my$bedmate1species1 = $ARGV[3];
my$bedmate2species1 = $ARGV[4];
my$bedmate1species2 = $ARGV[5];
my$bedmate2species2 = $ARGV[6];
#my$read_out = $ARGV[7]; 
#my$exon_out = $ARGV[8];
#my$gene_out = $ARGV[9];
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
    @elements = split(/\s+/,$_);
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

# initialize the database to store read names
#unlink "names" ;
system "touch $out_dir/names";
tie %readNames, "DB_File", "$out_dir/names", O_CREAT|O_RDWR, 0666, $DB_HASH or die "Cannot open file $out_dir/names: $!\n";

# initialize the database to store all read information (mate and species indicator)
#unlink "database" ;
system "touch $out_dir/database";
tie %readList, "DB_File", "$out_dir/database", O_CREAT|O_RDWR, 0666,$DB_HASH or die "Cannot open file $out_dir/database: $!\n";

# initialize the database to store file-specific reads to remove unique duplicates
#unlink "DB11" ;
system "touch $out_dir/DB11";
tie %readList11, "DB_File", "$out_dir/DB11", O_CREAT|O_RDWR, 0666,$DB_HASH or die "Cannot open file $out_dir/DB11: $!\n";

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
    @elements = split(/\s+/,$line);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene = $elements[0];
      $exon = $elements[1].",".$elements[2];
      $value = $gene."_".$exon; # CG12345_43567,98793
      $name = "1_1_".$read;
      #print "$name\t$value\n";
  
      if($readList{$name})
      {
        if((split("_",$readList{$name}))[0] eq $gene && (split("_",$readList{$name}))[1] eq $exon){next LINE;}
        else{delete $readList{$name}; next LINE;}
      }
      elsif($readList11{$name}){next LINE;}

      $readList11{$name} = 1;
      $readNames{$read} = 1;
      $readList{$name} = $value;
      #print "$name\t$readList{$name}\n";
      #printf("%s\t%s\t%s\n",$name,(split("_",$readList{$name}))[0],(split("_",$readList{$name}))[1]);
    }
  }
} 
close BED;
untie %readList11;

print "\nDone with 1st hash!\n";

#unlink "DB21" ;
system "touch $out_dir/DB21";
tie %readList21, "DB_File", "$out_dir/DB21",O_CREAT|O_RDWR, 0666, $DB_HASH or die "Cannot open file $out_dir/DB21: $!\n";

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
    @elements = split(/\s+/,$line);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene = $elements[0];
      $exon = $elements[1].",".$elements[2];
      $value = $gene."_".$exon; # CG12345_43567,98793
      $name = "2_1_".$read;
      #print "$name\t$value\n";

      if($readList{$name})
      {
        if((split("_",$readList{$name}))[0] eq $gene && (split("_",$readList{$name}))[1] eq $exon){next LINE;}
        else{delete $readList{$name}; next LINE;}
      }
      elsif($readList21{$name}){next LINE;}

      $readList21{$name} = 1;
      $readNames{$read} = 1;
      $readList{$name} = $value;
    }
  }
}
close BED;
untie %readList21;

print "\nDone with 2nd hash!\n";

#unlink "DB12" ;
system "touch $out_dir/DB12";
tie %readList12, "DB_File", "$out_dir/DB12", O_CREAT|O_RDWR, 0666,$DB_HASH or die "Cannot open file $out_dir/DB12: $!\n";

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
    @elements = split(/\s+/,$line);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene = $elements[0];
      $exon = $elements[1].",".$elements[2];
      $value = $gene."_".$exon; # CG12345_43567,98793
      $name = "1_2_".$read;
      #print "$name\t$value\n";

      if($readList{$name})
      {
        if((split("_",$readList{$name}))[0] eq $gene && (split("_",$readList{$name}))[1] eq $exon){next LINE;}
        else{delete $readList{$name}; next LINE;}
      }
      elsif($readList12{$name}){next LINE;}

      $readList12{$name} = 1;
      $readNames{$read} = 1;
      $readList{$name} = $value;
    }
  }
}
close BED;
untie %readList12;

print "\nDone with 3rd hash!\n";

#unlink "DB22" ;
system "touch $out_dir/DB22";
tie %readList22, "DB_File", "$out_dir/DB22",O_CREAT|O_RDWR, 0666, $DB_HASH or die "Cannot open file $out_dir/DB22: $!\n";

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
    @elements = split(/\s+/,$line);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene = $elements[0];
      $exon = $elements[1].",".$elements[2];
      $value = $gene."_".$exon; # CG12345_43567,98793
      $name = "2_2_".$read;
      #print "$name\t$value\n";

      if($readList{$name})
      {
        if((split("_",$readList{$name}))[0] eq $gene && (split("_",$readList{$name}))[1] eq $exon){next LINE;}
        else{delete $readList{$name}; next LINE;}
      }
      elsif($readList22{$name}){next LINE;}

      $readList22{$name} = 1;
      $readNames{$read} = 1;
      $readList{$name} = $value;
    }
  }
}
close BED;
untie %readList22;

print "\nDone with 4th and final hash!\n";

##### start counting #####
#there should be 16 possibilities for alignments of both mates to both species (4 files, 2 possibilities each, 2^4 = 16)
#these can be taken care of on a case by case basis

open(OUT1,">$out_dir\/$suffix\.mosaik.reads.txt") or die "Error writing to $out_dir\/$suffix\.mosaik.reads.txt\n";
print OUT1 "read\tgene1$species1\tgene1$species2\tgene2$species1\tgene2$species2\tcall"; 

foreach $read (keys %readNames)   
{   
  if($readList{'1_1_'.$read} && $readList{'1_2_'.$read}) # mate1 aligned to both genomes
  {
    if($readList{'2_1_'.$read} &&          #11
       $readList{'2_2_'.$read} && 
       (split("_",$readList{'1_1_'.$read}))[0] eq (split("_",$readList{'2_1_'.$read}))[0] && 
       (split("_",$readList{'1_2_'.$read}))[0] eq (split("_",$readList{'2_2_'.$read}))[0]) 
    {
      if((split("_",$readList{'1_1_'.$read}))[0] eq (split("_",$readList{'1_2_'.$read}))[0])
      {
        $geneHash{(split("_",$readList{'1_1_'.$read}))[0]}{b} += 1; 
        if((split("_",$readList{'1_1_'.$read}))[1] eq (split("_",$readList{'2_1_'.$read}))[1])
        {
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{b} += 1;
        } 
        else
        {
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{b} += 1; 
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'2_1_'.$read}))[1]}{b} += 1;
        }     
        printf OUT1 ("\n%s\t%s\t%s\t%s\t%s\tBoth",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'1_2_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0]);
      } 
      else{printf OUT1 ("\n%s\t%s\t%s\t%s\t%s\tError11",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'1_2_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0]);}
    } 
    elsif($readList{'2_1_'.$read}) #1
    {
      if((split("_",$readList{'1_1_'.$read}))[0] eq (split("_",$readList{'2_1_'.$read}))[0])
      {
        $geneHash{(split("_",$readList{'1_1_'.$read}))[0]}{$species1} += 1; 
        if((split("_",$readList{'1_1_'.$read}))[1] eq (split("_",$readList{'2_1_'.$read}))[1])
        {
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{$species1} += 1;
        } 
        else
        {
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{$species1} += 1; 
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'2_1_'.$read}))[1]}{$species1} += 1;
        } 
        printf OUT1 ("\n%s\t%s\t%s\t%s\t\*\t%s",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'1_2_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0],$species1);
      }
      else{printf OUT1 ("\n%s\t%s\t%s\t%s\t\*\tError1",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'1_2_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0])}
    }
    elsif($readList{'2_2_'.$read}) #6
    {
      if((split("_",$readList{'1_2_'.$read}))[0] eq (split("_",$readList{'2_2_'.$read}))[0])
      {
        $geneHash{(split("_",$readList{'1_2_'.$read}))[0]}{$species2} += 1; 
        if((split("_",$readList{'1_2_'.$read}))[1] eq (split("_",$readList{'2_2_'.$read}))[1])
        {
          $exonHash{(split("_",$readList{'1_2_'.$read}))[0]}{(split("_",$readList{'1_2_'.$read}))[1]}{$species2} += 1;
        } 
        else 
        {
          $exonHash{(split("_",$readList{'1_2_'.$read}))[0]}{(split("_",$readList{'1_2_'.$read}))[1]}{$species2} += 1; 
          $exonHash{(split("_",$readList{'1_2_'.$read}))[0]}{(split("_",$readList{'2_2_'.$read}))[1]}{$species2} += 1;
        } 
        printf OUT1 ("\n%s\t%s\t\*\t%s\t%s\t%s",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'1_2_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0],$species2);
      }
      else{printf OUT1 ("\n%s\t%s\t\*\t%s\t%s\tError6",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'1_2_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0]);}
    }
    else #13
    {
      if((split("_",$readList{'1_1_'.$read}))[0] eq (split("_",$readList{'1_2_'.$read}))[0] && 
         (split("_",$readList{'1_1_'.$read}))[1] eq (split("_",$readList{'1_2_'.$read}))[1])
      {
        $geneHash{(split("_",$readList{'1_1_'.$read}))[0]}{b} += 1; 
        $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{b} += 1;    
        printf OUT1 ("\n%s\t%s\t%s\t\*\t\*\tBoth",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'1_2_'.$read}))[0]);
      }       
      else{printf OUT1 ("\n%s\t%s\t%s\t\*\t\*\tError13",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'1_2_'.$read}))[0]);}
    } 
  }

  elsif($readList{'1_1_'.$read}) # mate1 aligned to mel but not sec
  {
    if($readList{'2_1_'.$read} && 
       $readList{'2_2_'.$read}) #2
    {
      if((split("_",$readList{'1_1_'.$read}))[0] eq (split("_",$readList{'2_1_'.$read}))[0])
      {
        $geneHash{(split("_",$readList{'1_1_'.$read}))[0]}{$species1} += 1; 
        if((split("_",$readList{'1_1_'.$read}))[1] eq (split("_",$readList{'2_1_'.$read}))[1])
        {
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{$species1} += 1;
        } 
        else
        {
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{$species1} += 1; 
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'2_1_'.$read}))[1]}{$species1} += 1;
        } 
        printf OUT1 ("\n%s\t%s\t\*\t%s\t%s\t%s",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0],$species1);
      }
      else{printf OUT1 ("\n%s\t%s\t\*\t%s\t%s\tError2",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0]);}
    } 
    elsif($readList{'2_1_'.$read}) #3
    {
      if((split("_",$readList{'1_1_'.$read}))[0] eq (split("_",$readList{'2_1_'.$read}))[0])
      {
        $geneHash{(split("_",$readList{'1_1_'.$read}))[0]}{$species1} += 1; 
        if((split("_",$readList{'1_1_'.$read}))[1] eq (split("_",$readList{'2_1_'.$read}))[1])
        {
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{$species1} += 1;
        } 
        else
        {
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{$species1} += 1; 
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'2_1_'.$read}))[1]}{$species1} += 1;
        } 
        printf OUT1 ("\n%s\t%s\t\*\t%s\t\*\t%s",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0],$species1);
      }
      else{printf OUT1 ("\n%s\t%s\t\*\t%s\t\*\tError3",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0]);}
    }
    elsif($readList{'2_2_'.$read}) #14
    {
      printf OUT1 ("\n%s\t%s\t\*\t\*\t%s\tTransplicing",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0]);
    }
    else #4
    {
      $geneHash{(split("_",$readList{'1_1_'.$read}))[0]}{$species1} += 1; 
      $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{$species1} += 1; 
      printf OUT1 ("\n%s\t%s\t\*\t\*\t\*\t%s",$read,(split("_",$readList{'1_1_'.$read}))[0],$species1);
    }       
  }

  elsif($readList{'1_2_'.$read}) # mate1 aligned to sec but not mel
  {
    if($readList{'2_1_'.$read} && 
       $readList{'2_2_'.$read}) #7
    {
      if((split("_",$readList{'1_2_'.$read}))[0] eq (split("_",$readList{'2_2_'.$read}))[0])
      {
        $geneHash{(split("_",$readList{'1_2_'.$read}))[0]}{$species2} += 1; 
        if((split("_",$readList{'1_2_'.$read}))[1] eq (split("_",$readList{'2_2_'.$read}))[1])
        {
          $exonHash{(split("_",$readList{'1_2_'.$read}))[0]}{(split("_",$readList{'1_2_'.$read}))[1]}{$species2} += 1;
        } 
        else
        {
          $exonHash{(split("_",$readList{'1_2_'.$read}))[0]}{(split("_",$readList{'1_2_'.$read}))[1]}{$species2} += 1; 
          $exonHash{(split("_",$readList{'1_2_'.$read}))[0]}{(split("_",$readList{'2_2_'.$read}))[1]}{$species2} += 1;
        }  
        printf OUT1 ("\n%s\t\*\t%s\t%s\t%s\t%s",$read,(split("_",$readList{'1_2_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0],$species2);
      }
      else{printf OUT1 ("\n%s\t\*\t%s\t%s\t%s\tError7",$read,(split("_",$readList{'1_2_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0]);}
    } 
    elsif($readList{'2_2_'.$read}) #8
    {
      if((split("_",$readList{'1_2_'.$read}))[0] eq (split("_",$readList{'2_2_'.$read}))[0])
      {
        $geneHash{(split("_",$readList{'1_2_'.$read}))[0]}{$species2} += 1; 
        if((split("_",$readList{'1_2_'.$read}))[1] eq (split("_",$readList{'2_2_'.$read}))[1])
        {
          $exonHash{(split("_",$readList{'1_2_'.$read}))[0]}{(split("_",$readList{'1_2_'.$read}))[1]}{$species2} += 1;
        } 
        else
        {
          $exonHash{(split("_",$readList{'1_2_'.$read}))[0]}{(split("_",$readList{'1_2_'.$read}))[1]}{$species2} += 1; 
          $exonHash{(split("_",$readList{'1_2_'.$read}))[0]}{(split("_",$readList{'2_2_'.$read}))[1]}{$species2} += 1;
        }  
        printf OUT1 ("\n%s\t\*\t%s\t\*\t%s\t%s",$read,(split("_",$readList{'1_2_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0],$species2);
      }
      else{printf OUT1 ("\n%s\t\*\t%s\t\*\t%s\tError8",$read,(split("_",$readList{'1_2_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0]);}
    }
    elsif($readList{'2_1_'.$read}) #15
    {
      printf OUT1 ("\n%s\t\*\t%s\t%s\t\*\tTransplicing",$read,(split("_",$readList{'1_2_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0]);
    } 
    else #9
    {
      $geneHash{(split("_",$readList{'1_2_'.$read}))[0]}{$species2} += 1; 
      $exonHash{(split("_",$readList{'1_2_'.$read}))[0]}{(split("_",$readList{'1_2_'.$read}))[1]}{$species2} += 1; 
      printf OUT1 ("\n%s\t\*\t%s\t\*\t\*\t%s",$read,(split("_",$readList{'1_2_'.$read}))[0],$species2);        
    }
  }

  else # mate1 aligned to nothing
  {
    if($readList{'2_1_'.$read} && 
       $readList{'2_2_'.$read}) #12
    {
      if((split("_",$readList{'2_1_'.$read}))[0] eq (split("_",$readList{'2_2_'.$read}))[0] && 
         (split("_",$readList{'2_1_'.$read}))[1] eq (split("_",$readList{'2_2_'.$read}))[1])
      {
        $geneHash{(split("_",$readList{'2_1_'.$read}))[0]}{b} += 1; 
        $exonHash{(split("_",$readList{'2_1_'.$read}))[0]}{(split("_",$readList{'2_1_'.$read}))[1]}{b} += 1;   
        printf OUT1 ("\n%s\t\*\t\*\t%s\t%s\tBoth",$read,(split("_",$readList{'2_1_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0]);
      }
      else{printf OUT1 ("\n%s\t\*\t\*\t%s\t%s\tError12",$read,(split("_",$readList{'2_1_'.$read}))[0],(split("_",$readList{'2_2_'.$read}))[0]);}
    }
    elsif($readList{'2_1_'.$read}) #5
    {
      $geneHash{(split("_",$readList{'2_1_'.$read}))[0]}{$species1} += 1; 
      $exonHash{(split("_",$readList{'2_1_'.$read}))[0]}{(split("_",$readList{'2_1_'.$read}))[1]}{$species1} += 1;   
      printf OUT1 ("\n%s\t\*\t\*\t%s\t\*\t%s",$read,(split("_",$readList{'2_1_'.$read}))[0],$species1);              
    } 
    elsif($readList{'2_2_'.$read}) #10
    {
      $geneHash{(split("_",$readList{'2_2_'.$read}))[0]}{$species2} += 1; 
      $exonHash{(split("_",$readList{'2_2_'.$read}))[0]}{(split("_",$readList{'2_2_'.$read}))[1]}{$species2} += 1;
      printf OUT1 ("\n%s\t\*\t\*\t%s\t\*\t%s",$read,(split("_",$readList{'2_2_'.$read}))[0],$species2);
    }
    else #16
    {
      printf OUT1 "\n$read\t*\t*\t*\t*\tNA";
    }       
  }
}
close OUT1;
untie %readNames;
untie %readList;

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
