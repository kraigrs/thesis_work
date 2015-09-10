#!/usr/bin/perl

###################################################################################
#
# 04/21/2010
#
# classify4.pl
#
# Purpose: (adapted from classify3.pl) 
# 
# Input: two separate alignment files (1 for each species)
#
# Output: a list of genes and their respective match and mismatch counts, along with a 
#         summary of all the reads 
#
# E.g. perl classify4.pl Dmel Dsec Hyb.Dmel_genes.dm3.map Hyb.Dsec_genes.droSec1.map.lifted Hyb.reads.txt Hyb.exons.txt Hyb.genes.txt &
#
###################################################################################

use strict;
use warnings;

my$line1; my$line2;
my$gene1mel; my$gene1sec; my$gene2mel; my$gene2sec;
my$exon1mel; my$exon1sec; my$exon2mel; my$exon2sec;
my$read1mel; my$read1sec; my$read2mel; my$read2sec;
my@elements;
my%geneHash;
my%exonHash;

my$species1 = $ARGV[0]; # Dmel
my$species2 = $ARGV[1]; # Dsec
my$BEDfile1 = $ARGV[2];
my$BEDfile2 = $ARGV[3];
my$read_out = $ARGV[4]; 

open(OUT1,">$read_out") or die "Error writing to $read_out\n";
print OUT1 "read\tgene1$species1\tgene1$species2\tgene2$species1\tgene2$species2\tcall";

#print "Printed first line\n";

open(BED1,"$BEDfile1") or die "Can't open $BEDfile1 for reading!";
open(BED2,"$BEDfile2") or die "Can't open $BEDfile2 for reading!";

LINE:while(<BED1>) 
{
  $line1 = $_;
  $line2 = <BED2>; 

  chomp $line1;
  chomp $line2;  
  #print "new line\n";
  if($line1 =~ /^\@.+/){next LINE;}

  #find mate 1
  elsif($line1 =~ /HWI.+\/1/)
  { 
    @elements = split(/\s/,$line1);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read1mel = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene1mel = $elements[0];
      $exon1mel = $elements[1] . "," . $elements[2];
    }    
    else
    {
      $gene1mel = "*";
      $exon1mel = "*";
    }

    @elements = split(/\s/,$line2);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read1sec = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene1sec = $elements[0];
      $exon1sec = $elements[1].",".$elements[2];
    }    
    else
    {
      $gene1sec = "*";
      $exon1sec = "*";
    }
    next LINE;
  }

  #find mate 2
  elsif($line1 =~ /HWI.+\/2/)
  { 
    @elements = split(/\s/,$line1);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read2mel = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene2mel = $elements[0];
      $exon2mel = $elements[1] . "," . $elements[2];
    }    
    else
    {
      $gene2mel = "*";
      $exon2mel = "*";
    }

    @elements = split(/\s/,$line2);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read2sec = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene2sec = $elements[0];
      $exon2sec = $elements[1].",".$elements[2];
    }    
    else
    {
      $gene2sec = "*";
      $exon2sec = "*";
    }

    if($read2mel ne $read2sec){die "Files do not line up!!!\n";}  
   
    unless($read1mel ne $read2mel && $read1sec ne $read2sec)
    {
      if($gene1mel ne "*" && $gene1sec ne "*") # mate1 aligned to both genomes
      {
        if($gene2mel ne "*" && $gene2sec ne "*" && $gene1mel eq $gene2mel && $gene1sec eq $gene2sec) #11
        {
          if($gene1mel eq $gene1sec)
	  {
            $geneHash{$gene1mel}{b} += 1; 
            if($exon1mel eq $exon2mel){$exonHash{$gene1mel}{$exon1mel}{b} += 1;} 
            else{$exonHash{$gene1mel}{$exon1mel}{b} += 1; $exonHash{$gene1mel}{$exon2mel}{b} += 1;}     
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tBoth";
            next LINE;
          } 
          else
          {
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tError11";
            next LINE;
          }
        } 
        elsif($gene2mel eq "*" && $gene2sec eq "*") #13
        {
          if($gene1mel eq $gene1sec)
	  {
            $geneHash{$gene1mel}{b} += 1; 
            $exonHash{$gene1mel}{$exon1mel}{b} += 1;    
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tBoth";
            next LINE;
          }     
          else
          {
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tError13";
            next LINE;
          }
        } 
        elsif($gene2mel ne "*" && $gene2sec eq "*" && $gene1mel eq $gene2mel) #1
        {
          if($gene1mel =~ /^C\w{1}\d+/)
	  {
            $geneHash{$gene1mel}{$species1} += 1; 
            if($exon1mel eq $exon2mel){$exonHash{$gene1mel}{$exon1mel}{$species1} += 1;} 
            else{$exonHash{$gene1mel}{$exon1mel}{$species1} += 1; $exonHash{$gene1mel}{$exon2mel}{$species1} += 1;} 
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\t$species1";
            next LINE;
          }
          else
          {
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tError1";
            next LINE;
          }
        }
        elsif($gene2mel eq "*" && $gene2sec ne "*" && $gene1sec eq $gene2sec) #6
        {
          $geneHash{$gene1sec}{$species2} += 1; 
          if($exon1sec eq $exon2sec){$exonHash{$gene1sec}{$exon1sec}{$species2} += 1;} 
          else{$exonHash{$gene1sec}{$exon1sec}{$species2} += 1; $exonHash{$gene1sec}{$exon2sec}{$species2} += 1;} 
          print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\t$species2";
          next LINE;
        }
        else
        {
          print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tInvestigate";
          next LINE;
        } 
      }

      elsif($gene1mel ne "*" && $gene1sec eq "*") # mate1 aligned to mel but not sec
      {
        if($gene2mel ne "*" && $gene1mel eq $gene2mel) #2 & 3
        {
          if($gene1mel =~ /^C\w{1}\d+/)
	  {
            $geneHash{$gene1mel}{$species1} += 1; 
            if($exon1mel eq $exon2mel){$exonHash{$gene1mel}{$exon1mel}{$species1} += 1;} 
            else{$exonHash{$gene1mel}{$exon1mel}{$species1} += 1; $exonHash{$gene1mel}{$exon2mel}{$species1} += 1;} 
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\t$species1";
            next LINE;
          }
          else
          {
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tError2and3";
            next LINE;
          }
        } 
        elsif($gene2mel eq "*" && $gene2sec eq "*") #4
        {
          if($gene1mel =~ /^C\w{1}\d+/)
	  {
            $geneHash{$gene1mel}{$species1} += 1; 
            $exonHash{$gene1mel}{$exon1mel}{$species1} += 1; 
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\t$species1";
            next LINE;
          } 
          else
          {
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tError4";
            next LINE;
          }
        }
        elsif($gene2mel eq "*" && $gene2sec ne "*") #14
        {
          print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tTransplicing";
          next LINE;
        }
        else
        {
          print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tInvestigate";
          next LINE;
        } 
      }

      elsif($gene1mel eq "*" && $gene1sec ne "*") # mate1 aligned to sec but not mel
      {
        if($gene2sec ne "*" && $gene1sec eq $gene2sec) #7 & 8
        {
          $geneHash{$gene1sec}{$species2} += 1; 
          if($exon1sec eq $exon2sec){$exonHash{$gene1sec}{$exon1sec}{$species2} += 1;} 
          else{$exonHash{$gene1sec}{$exon1sec}{$species2} += 1; $exonHash{$gene1sec}{$exon2sec}{$species2} += 1;}  
          print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\t$species2";
          next LINE;
        } 
        elsif($gene2mel eq "*" && $gene2sec eq "*") #9
        {
          $geneHash{$gene1sec}{$species2} += 1; 
          $exonHash{$gene1sec}{$exon1sec}{$species2} += 1; 
          print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\t$species2";
          next LINE;
        }
        elsif($gene2mel ne "*" && $gene2sec eq "*") #15
        {
          print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tTransplicing";
          next LINE;
        } 
        else
        {
          print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tInvestigate";
          next LINE;
        }
      }

      elsif($gene1mel eq "*" && $gene1sec eq "*") # mate1 aligned to nothing
      {
        if($gene2mel ne "*" && $gene2sec eq "*") #5
        {
          if($gene2mel =~ /^C\w{1}\d+/)
	  {
            $geneHash{$gene2mel}{$species1} += 1; 
            $exonHash{$gene2mel}{$exon2mel}{$species1} += 1;   
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\t$species1";
            next LINE;
          }     
          else
          {
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tError5";
            next LINE;
          }
        } 
        elsif($gene2mel eq "*" && $gene2sec ne "*") #10
        {
          $geneHash{$gene2sec}{$species2} += 1; 
          $exonHash{$gene2sec}{$exon2sec}{$species2} += 1;
          print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\t$species2";
          next LINE;
        }
        elsif($gene2mel ne "*" && $gene2sec ne "*") #12
        {
          if($gene2mel eq $gene2sec)
	  {
            $geneHash{$gene2mel}{b} += 1; 
            $exonHash{$gene2mel}{$exon2mel}{b} += 1;   
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tBoth";
            next LINE;
	  }
          else
          {
            print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tError10";
            next LINE;
          }
        } 
        elsif($gene2mel eq "*" && $gene2sec eq "*") #16
        {
          print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tNA";
          next LINE;
        }
        else
        {
          print OUT1 "\n$read1mel\t$gene1mel\t$gene1sec\t$gene2mel\t$gene2sec\tInvestigate";
          next LINE;
        }
      }
    }       
  }
}
close BED1;
close BED2;
close OUT1;


my$exon_out = $ARGV[5];
open(OUT2,"> $exon_out") or die "Error writing to $exon_out\n";
print OUT2 "gene\texon\t$species1\t$species2\tBoth\tAdj$species1\tAdj$species2";

my$madj = 0;
my$sadj = 0;

foreach my$k1 (sort keys %exonHash)
{
  foreach my$k2 (sort keys %{$exonHash{$k1}})
  {
    if($exonHash{$k1}{$k2}{$species1} && $exonHash{$k1}{$k2}{$species2} && $exonHash{$k1}{$k2}{b})
    {
      $madj = $exonHash{$k1}{$k2}{$species1} + ($exonHash{$k1}{$k2}{$species1} / ($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2})) * $exonHash{$k1}{$k2}{b};  
      $sadj = $exonHash{$k1}{$k2}{$species2} + ($exonHash{$k1}{$k2}{$species2} / ($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2})) * $exonHash{$k1}{$k2}{b};

      print OUT2 "\n$k1\t$k2\t$exonHash{$k1}{$k2}{$species1}\t$exonHash{$k1}{$k2}{$species2}\t$exonHash{$k1}{$k2}{b}\t$madj\t$sadj";
    }
    elsif($exonHash{$k1}{$k2}{$species1} && $exonHash{$k1}{$k2}{$species2})
    {
      print OUT2 "\n$k1\t$k2\t$exonHash{$k1}{$k2}{$species1}\t$exonHash{$k1}{$k2}{$species2}\t0\t$exonHash{$k1}{$k2}{$species1}\t$exonHash{$k1}{$k2}{$species2}";
    }
    elsif($exonHash{$k1}{$k2}{$species1} && $exonHash{$k1}{$k2}{b})
    {
      $madj = $exonHash{$k1}{$k2}{$species1} + 0.5 * $exonHash{$k1}{$k2}{b};  

      print OUT2 "\n$k1\t$k2\t$exonHash{$k1}{$k2}{$species1}\t0\t$exonHash{$k1}{$k2}{b}\t$madj\t0";
    }
    elsif($exonHash{$k1}{$k2}{$species2} && $exonHash{$k1}{$k2}{b})
    {
      $sadj = $exonHash{$k1}{$k2}{$species2} + 0.5 * $exonHash{$k1}{$k2}{b};

      print OUT2 "\n$k1\t$k2\t0\t$exonHash{$k1}{$k2}{$species2}\t$exonHash{$k1}{$k2}{b}\t0\t$sadj";
    }
    elsif($exonHash{$k1}{$k2}{$species1})
    {
      print OUT2 "\n$k1\t$k2\t$exonHash{$k1}{$k2}{$species1}\t0\t0\t$exonHash{$k1}{$k2}{$species1}\t0";
    }
    elsif($exonHash{$k1}{$k2}{$species2})
    {
      print OUT2 "\n$k1\t$k2\t0\t$exonHash{$k1}{$k2}{$species2}\t0\t0\t$exonHash{$k1}{$k2}{$species2}";
    }
    elsif($exonHash{$k1}{$k2}{b})
    {
      print OUT2 "\n$k1\t$k2\t0\t0\t$exonHash{$k1}{$k2}{b}\t0\t0";
    }
    else
    {
      print OUT2 "\n$k1\t$k2\t0\t0\t0";
    }
    $madj = 0;
    $sadj = 0;
  }
}
close OUT2;

my$gene_out = $ARGV[6];
open(OUT3,"> $gene_out") or die "Error writing to $gene_out\n";
print OUT3 "gene\t$species1\t$species2\tBoth\tAdj$species1\tAdj$species2";

$madj = 0;
$sadj = 0;

foreach(sort keys %geneHash)
{
  if($geneHash{$_}{$species1} && $geneHash{$_}{$species2} && $geneHash{$_}{b})
  {
    $madj = $geneHash{$_}{$species1} + ($geneHash{$_}{$species1} / ($geneHash{$_}{$species1} + $geneHash{$_}{$species2})) * $geneHash{$_}{b};  
    $sadj = $geneHash{$_}{$species2} + ($geneHash{$_}{$species2} / ($geneHash{$_}{$species1} + $geneHash{$_}{$species2})) * $geneHash{$_}{b};

    print OUT3 "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}\t$geneHash{$_}{b}\t$madj\t$sadj";
  }
  elsif($geneHash{$_}{$species1} && $geneHash{$_}{$species2})
  {
    print OUT3 "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}\t0\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}";
  }
  elsif($geneHash{$_}{$species1} && $geneHash{$_}{b})
  {
    $madj = $geneHash{$_}{$species1} + 0.5 * $geneHash{$_}{b};  
  
    print OUT3 "\n$_\t$geneHash{$_}{$species1}\t0\t$geneHash{$_}{b}\t$madj\t0";
  }
  elsif($geneHash{$_}{$species2} && $geneHash{$_}{b})
  {
    $sadj = $geneHash{$_}{$species2} + 0.5 * $geneHash{$_}{b};

    print OUT3 "\n$_\t0\t$geneHash{$_}{$species2}\t$geneHash{$_}{b}\t0\t$sadj";
  }
  elsif($geneHash{$_}{$species1})
  {
    print OUT3 "\n$_\t$geneHash{$_}{$species1}\t0\t0\t$geneHash{$_}{$species1}\t0";
  }
  elsif($geneHash{$_}{$species2})
  {
    print OUT3 "\n$_\t0\t$geneHash{$_}{$species2}\t0\t0\t$geneHash{$_}{$species2}";
  }
  elsif($geneHash{$_}{b})
  {
    print OUT3 "\n$_\t0\t0\t$geneHash{$_}{b}\t0\t0";
  }
  else
  {
    print OUT3 "\n$_\t0\t0\t0";
  }
  $madj = 0;
  $sadj = 0;
}
close OUT3;
