#!/usr/bin/perl

###################################################################################
#
# 03/21/2011
#
# classify.0.6.pl
#
# Purpose: combine 2 .bed files (lifted) to elucidate allele-specific expression. the alignment files are
#          from individual species aligned to their respective genomes
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
# E.g. perl classify.0.6.pl <const> <species> <bedmate1species> <bedmate2species> <readout> <exonout> <geneout>
#     
#      <const>            ==> list of constitutive exons
#      <species>          ==> indicate name of species aligned to (i.e. Dmel, Dsec, etc.)
#      <bedmate1species>  ==> file with mate 1 alignments to species
#      <bedmate2species>  ==> file with mate 2 alignments to species
#      <readout>          ==> file to send reads for output
#      <exonout>          ==> file to send exon counts for output
#      <geneout>          ==> file to send gene counts for output
#
# perl classify.0.6.pl ../McManus/const_exons.dm3.bed zhr ../mel_sim_data/mRNA-Seq/zhr/zhr.mate1.zhr.mosaik.bed.lifted.converted ../mel_sim_data/mRNA-Seq/zhr/zhr.mate2.zhr.mosaik.bed.lifted.converted ../mel_sim_data/mRNA-Seq/zhr/zhr.mRNA.mosaik.reads.txt ../mel_sim_data/mRNA-Seq/zhr/zhr.mRNA.mosaik.exons.txt ../mel_sim_data/mRNA-Seq/zhr/zhr.mRNA.mosaik.genes.txt &
#
# perl classify.0.6.pl ../McManus/constitutive_regions.bed zhr ../mel_sim_data/mRNA-Seq/zhr/zhr.mate1.zhr.mosaik.bed.lifted.converted ../mel_sim_data/mRNA-Seq/zhr/zhr.mate2.zhr.mosaik.bed.lifted.converted ./zhr.mRNA.mosaik.reads.txt ./zhr.mRNA.mosaik.exons.txt ./zhr.mRNA.mosaik.genes.txt &
#
# perl classify.0.6.pl ../McManus/const_exons.dm3.bed tsimbazaza ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate1.tsimbazaza.mosaik.bed.lifted.converted ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mate2.tsimbazaza.mosaik.bed.lifted.converted ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mRNA.mosaik.reads.txt ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mRNA.mosaik.exons.txt ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mRNA.mosaik.genes.txt &
#
###################################################################################

use strict;
use warnings;
use DB_File;
use Fcntl;

main();
sub main
{

my$const = $ARGV[0];
my$species = $ARGV[1]; # Dmel
my$bedmate1species = $ARGV[2];
my$bedmate2species = $ARGV[3];
my$read_out = $ARGV[4]; 
my$exon_out = $ARGV[5];
my$gene_out = $ARGV[6];

our($line,$value,$gene,$exon,$read,$name);
my@elements;
our(%geneHash,%exonHash,%readNames,%readList,%readList11,%readList21);

# initialize the database to store read names
#unlink "names";
system "touch names";
tie %readNames,"DB_File","names",$DB_HASH or die "Cannot open file 'names': $!\n";

# initialize the database to store all read information (mate and species indicator)
#unlink "database";
system "touch database";
tie %readList,"DB_File","database",$DB_HASH or die "Cannot open file 'database': $!\n";

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

    $geneHash{$gene}{$species} = 0; 
    $exonHash{$gene}{$exon}{$species} = 0;
  }
}
close CONST;
print "\nDone with genes and exons hashes!\n";

# initialize the database to store file-specific reads to remove unique duplicates
#unlink "DB11";
system "touch DB11";
tie %readList11,"DB_File","DB11",$DB_HASH or die "Cannot open file 'DB11': $!\n";

# mate 1
open(BED,"$bedmate1species") or die "Can't open $bedmate1species for reading!";
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

#unlink "DB21";
system "touch DB21";
tie %readList21,"DB_File","DB21",$DB_HASH or die "Cannot open file 'DB21': $!\n";

# mate 2, species 1
open(BED,"$bedmate2species") or die "Can't open $bedmate2species for reading!";
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

print "\nDone with hashes!\n";

##### start counting #####

open(OUT1,">$read_out") or die "Error writing to $read_out\n";
print OUT1 "read\tgene1$species\tgene2$species\tcall"; 

foreach $read (keys %readNames)   
{   

  if($readList{'1_1_'.$read}) # mate1 aligned to genome
  {
    if($readList{'2_1_'.$read}) # mate2 also aligned to genome
    {
      if((split("_",$readList{'1_1_'.$read}))[0] eq (split("_",$readList{'2_1_'.$read}))[0]) # mate1 and mate2 align to same gene!
      {
        $geneHash{(split("_",$readList{'1_1_'.$read}))[0]}{$species} += 1; 
        if((split("_",$readList{'1_1_'.$read}))[1] eq (split("_",$readList{'2_1_'.$read}))[1])
        {
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{$species} += 1;
        } 
        else
        {
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{$species} += 1; 
          $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'2_1_'.$read}))[1]}{$species} += 1;
        } 
        printf OUT1 ("\n%s\t%s\t%s\t%s",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0],$species);
      }
      else{printf OUT1 ("\n%s\t%s\t%s\tError",$read,(split("_",$readList{'1_1_'.$read}))[0],(split("_",$readList{'2_1_'.$read}))[0]);}
    } 
    else # mate2 did not align to genome
    {
      $geneHash{(split("_",$readList{'1_1_'.$read}))[0]}{$species} += 1; 
      $exonHash{(split("_",$readList{'1_1_'.$read}))[0]}{(split("_",$readList{'1_1_'.$read}))[1]}{$species} += 1; 
      printf OUT1 ("\n%s\t%s\t\*\t%s",$read,(split("_",$readList{'1_1_'.$read}))[0],$species);
    }       
  }

  else # mate1 aligned to nothing
  {
    if($readList{'2_1_'.$read}) # ... but mate2 aligned!
    {
      $geneHash{(split("_",$readList{'2_1_'.$read}))[0]}{$species} += 1; 
      $exonHash{(split("_",$readList{'2_1_'.$read}))[0]}{(split("_",$readList{'2_1_'.$read}))[1]}{$species} += 1;   
      printf OUT1 ("\n%s\t\*\t%s\t%s",$read,(split("_",$readList{'2_1_'.$read}))[0],$species);              
    }
    else # nothing aligned
    {
      print OUT1 "\n$read\t*\t*\tNA";
    }       
  }
}
close OUT1;
untie %readNames;
untie %readList;

open(OUT2,"> $exon_out") or die "Error writing to $exon_out\n";
print OUT2 "gene_exon\t$species";

foreach my$k1 (keys %exonHash)
{
  foreach my$k2 (keys %{$exonHash{$k1}})
  {
    if($exonHash{$k1}{$k2}{$species})
    {
      print OUT2 "\n$k1\_$k2\t$exonHash{$k1}{$k2}{$species}";
    }
  }
}
close OUT2;

open(OUT3,"> $gene_out") or die "Error writing to $gene_out\n";
print OUT3 "gene\t$species";

foreach(keys %geneHash)
{
  if($geneHash{$_}{$species})
  {
    print OUT3 "\n$_\t$geneHash{$_}{$species}";
  }
}
close OUT3;
}
