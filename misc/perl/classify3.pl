#!/usr/bin/perl

###################################################################################
#
# 03/03/2009
#
# classify3.pl
#
# Purpose: (adapted from classify2.pl) classify reads (single and paired-end) as matching
#          or mismatching aligned separately to reference genomes. The difference here is not using 
#          SNP information, but rather simply aligning to a particular species in the reference
#          i.e. the reference contains D. mel. and D. sec. references, whichever a read aligns to,
#          it belongs to that species. this is complicated by the fact that we are aligning each
#          read independently to two references, so there is the potential to align to multiple 
#          genes and/or species for any given mate paired-end read
# 
# Input: two separate alignment files (mates 1 and mates 2) 
#
# Output: a list of genes and their respective match and mismatch counts, along with a 
#         summary of all the reads 
#
# E.g. perl classify3.pl names_mel-sec.txt Dmel Dsec Hyb.mate1.Dmel-exon.sam Hyb.mate1.Dsec-exon.sam Hyb.mate2.Dmel-exon.sam Hyb.mate2.Dsec-exon.sam Hyb.reads.txt Hyb.exons.txt Hyb.genes.txt &
#
###################################################################################

use strict;
#use warnings;

my@elements;
my%sec2mel;
my%mel2CG;

my$annotation = $ARGV[0];
open(NAMES,"$annotation") or die "Can't open $annotation for reading!";
while(<NAMES>)
{
  @elements = split(/\s+/,$_); 

  $sec2mel{$elements[0]} = $elements[2];
  $mel2CG{$elements[1]} = $elements[2];  
}
close NAMES;

my$line1mel; my$line1sec; my$line2mel; my$line2sec;
my$read1mel; my$read1sec; my$read2mel; my$read2sec;
my$align1mel; my$align1sec; my$align2mel; my$align2sec;
my$gene1mel; my$gene1sec; my$gene2mel; my$gene2sec;
my$spec1mel; my$spec1sec; my$spec2mel; my$spec2sec;
my$exon1mel; my$exon1sec; my$exon2mel; my$exon2sec;
my$seq1mel; my$seq1sec; my$seq2mel; my$seq2sec;
my@elements1mel; my@elements1sec; my@elements2mel; my@elements2sec;
my%geneHash;
my%exonHash;

my$species1 = $ARGV[1]; # Dmel
my$species2 = $ARGV[2]; # Dsec
my$SAMfile1mel = $ARGV[3];
my$SAMfile1sec = $ARGV[4];
my$SAMfile2mel = $ARGV[5];
my$SAMfile2sec = $ARGV[6];
my$read_out = $ARGV[7]; 

open(OUT,">> $read_out") or die "Error writing to $read_out\n";
print OUT "read\tseq1\tseq2\tmate1$species1\tmate1$species2\tmate2$species1\tmate2$species2\tcall";
close OUT;

open(SAM1,"$SAMfile1mel") or die "Can't open $SAMfile1mel for reading!";
open(SAM2,"$SAMfile1sec") or die "Can't open $SAMfile1sec for reading!";
open(SAM3,"$SAMfile2mel") or die "Can't open $SAMfile2mel for reading!";
open(SAM4,"$SAMfile2sec") or die "Can't open $SAMfile2sec for reading!";
while (<SAM1>) 
{
  chomp;
  $line1mel = $_;
  $line1sec = <SAM2>; 
  $line2mel = <SAM3>;
  $line2sec = <SAM4>;
  
  if($line1mel =~ /^\@.+/){next;}
  elsif($line1mel =~ /^HWI.+\/1\s+/)
  { 
    @elements1mel = split(/\s/,$line1mel);
    if($elements1mel[0] =~ /^(HWI.+)\/\d{1}/){$read1mel = $1;}
    $align1mel = $elements1mel[2];
    if($align1mel =~ /^(.+)\:(.+)$/)
    {
      $gene1mel = $1; $exon1mel = $2;
      #print "\nmate1\t$gene1mel\t$exon1mel";
    }
    else{$gene1mel = "*"; $exon1mel = "*";}
    $seq1mel = $elements1mel[9];
    #print "$seq1mel\t";

    @elements1sec = split(/\s/,$line1sec);
    if($elements1sec[0] =~ /^(HWI.+)\/\d{1}/){$read1sec = $1;}
    $align1sec = $elements1sec[2];
    if($align1sec =~ /^(.+)\:(.+)$/)
    {
      $gene1sec = $1; $exon1sec = $2;
      #print "\nmate1\t$gene1sec\t$exon1sec";
    }
    else{$gene1sec = "*"; $exon1sec = "*";}
    $seq1sec = $elements1sec[9];
    #print "$seq1sec\t";
  
    @elements2mel = split(/\s/,$line2mel);
    if($elements2mel[0] =~ /^(HWI.+)\/\d{1}/){$read2mel = $1;}
    $align2mel = $elements2mel[2];
    if($align2mel =~ /^(.+)\:(.+)$/)
    {
      $gene2mel = $1; $exon2mel = $2;
      #print "\nmate2\t$gene2mel\t$exon2mel";
    }
    else{$gene2mel = "*"; $exon2mel = "*";}
    $seq2mel = $elements2mel[9];
    #print "$seq2mel\n";

    @elements2sec = split(/\s/,$line2sec);
    if($elements2sec[0] =~ /^(HWI.+)\/\d{1}/){$read2sec = $1;}
    $align2sec = $elements2sec[2];
    if($align2sec =~ /^(.+)\:(.+)$/)
    {
      $gene2sec = $1; $exon2sec = $2;
      #print "\nmate2\t$gene2sec\t$exon2sec";
    }
    else{$gene2sec = "*"; $exon2sec = "*";}
    $seq2sec = $elements2sec[9];
    #print "$seq2sec\n";
   
    unless($read1mel ne $read2mel && $read1sec ne $read2sec)
    {
      if($gene1mel ne "*" && $gene1sec ne "*") # mate1 aligned to both genomes
      {
        if($gene2mel ne "*" && $gene2sec ne "*" && $gene1mel eq $gene2mel && $gene1sec eq $gene2sec && $sec2mel{$gene1sec}) #11
        {
          if($gene1mel eq $sec2mel{$gene1sec})
	  {
            $geneHash{$gene1mel}{b} += 1; 
            if($exon1mel == $exon2mel){$exonHash{$gene1mel}{$exon1mel}{b} += 1;} 
            else{$exonHash{$gene1mel}{$exon1mel}{b} += 1; $exonHash{$gene1mel}{$exon2mel}{b} += 1;}     
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tBoth";
            close OUT; 
          } 
          elsif($mel2CG{$gene1mel} eq $sec2mel{$gene1sec})
	  {
            $geneHash{$mel2CG{$gene1mel}}{b} += 1; 
            if($exon1mel == $exon2mel){$exonHash{$mel2CG{$gene1mel}}{$exon1mel}{b} += 1;} 
            else{$exonHash{$mel2CG{$gene1mel}}{$exon1mel}{b} += 1; $exonHash{$mel2CG{$gene1mel}}{$exon2mel}{b} += 1;}     
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tBoth";
            close OUT;
          }
          else
          {
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tError11";
            close OUT;       
          }
        } 
        elsif($gene2mel eq "*" && $gene2sec eq "*" && $sec2mel{$gene1sec}) #13
        {
          if($gene1mel eq $sec2mel{$gene1sec})
	  {
            $geneHash{$gene1mel}{b} += 1; 
            $exonHash{$gene1mel}{$exon1mel}{b} += 1;    
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tBoth";
            close OUT;
          }    
          elsif($mel2CG{$gene1mel} eq $sec2mel{$gene1sec})
	  {
            $geneHash{$mel2CG{$gene1mel}}{b} += 1; 
            $exonHash{$mel2CG{$gene1mel}}{$exon1mel}{b} += 1;    
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tBoth";
            close OUT; 
          } 
          else
          {
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tError13";
            close OUT;       
          }
        } 
        elsif($gene2mel ne "*" && $gene2sec eq "*" && $gene1mel eq $gene2mel) #1
        {
          if($gene1mel =~ /^C\w{1}\d+/)
	  {
            $geneHash{$gene1mel}{$species1} += 1; 
            if($exon1mel == $exon2mel){$exonHash{$gene1mel}{$exon1mel}{$species1} += 1;} 
            else{$exonHash{$gene1mel}{$exon1mel}{$species1} += 1; $exonHash{$gene1mel}{$exon2mel}{$species1} += 1;} 
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\t$species1";
            close OUT;
          }
          elsif($mel2CG{$gene1mel})
	  {
            $geneHash{$mel2CG{$gene1mel}}{$species1} += 1; 
            if($exon1mel == $exon2mel){$exonHash{$mel2CG{$gene1mel}}{$exon1mel}{$species1} += 1;} 
            else{$exonHash{$mel2CG{$gene1mel}}{$exon1mel}{$species1} += 1; $exonHash{$mel2CG{$gene1mel}}{$exon2mel}{$species1} += 1;} 
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\t$species1";
            close OUT;
            #print "$mel2CG{$gene1mel}\n";
          }
          else
          {
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tError1";
            close OUT;       
          }
        }
        elsif($gene2mel eq "*" && $gene2sec ne "*" && $gene1sec eq $gene2sec && $sec2mel{$gene1sec}) #6
        {
          $geneHash{$sec2mel{$gene1sec}}{$species2} += 1; 
          if($exon1sec == $exon2sec){$exonHash{$sec2mel{$gene1sec}}{$exon1sec}{$species2} += 1;} 
          else{$exonHash{$sec2mel{$gene1sec}}{$exon1sec}{$species2} += 1; $exonHash{$sec2mel{$gene1sec}}{$exon2sec}{$species2} += 1;} 
          open(OUT,">> $read_out") or die "Error writing to $read_out\n";
          print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\t$species2";
          close OUT;       
        }
        else
        {
          open(OUT,">> $read_out") or die "Error writing to $read_out\n";
          print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tInvestigate";
          close OUT;       
        } 
      }

      elsif($gene1mel ne "*" && $gene1sec eq "*") # mate1 aligned to mel but not sec
      {
        if($gene2mel ne "*" && $gene1mel eq $gene2mel) #2 & 3
        {
          if($gene1mel =~ /^C\w{1}\d+/)
	  {
            $geneHash{$gene1mel}{$species1} += 1; 
            if($exon1mel == $exon2mel){$exonHash{$gene1mel}{$exon1mel}{$species1} += 1;} 
            else{$exonHash{$gene1mel}{$exon1mel}{$species1} += 1; $exonHash{$gene1mel}{$exon2mel}{$species1} += 1;} 
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\t$species1";
            close OUT;         
          }
          elsif($mel2CG{$gene1mel})
	  {
            $geneHash{$mel2CG{$gene1mel}}{$species1} += 1; 
            if($exon1mel == $exon2mel){$exonHash{$mel2CG{$gene1mel}}{$exon1mel}{$species1} += 1;} 
            else{$exonHash{$mel2CG{$gene1mel}}{$exon1mel}{$species1} += 1; $exonHash{$mel2CG{$gene1mel}}{$exon2mel}{$species1} += 1;} 
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\t$species1";
            close OUT;         
          }
          else
          {
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tError2and3";
            close OUT;       
          }
        } 
        elsif($gene2mel eq "*" && $gene2sec eq "*") #4
        {
          if($gene1mel =~ /^C\w{1}\d+/)
	  {
            $geneHash{$gene1mel}{$species1} += 1; 
            $exonHash{$gene1mel}{$exon1mel}{$species1} += 1; 
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\t$species1";
            close OUT;
          } 
          elsif($mel2CG{$gene1mel})
	  {
            $geneHash{$mel2CG{$gene1mel}}{$species1} += 1; 
            $exonHash{$mel2CG{$gene1mel}}{$exon1mel}{$species1} += 1; 
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\t$species1";
            close OUT;
          }
          else
          {
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tError4";
            close OUT;       
          }
        }
        elsif($gene2mel eq "*" && $gene2sec ne "*") #14
        {
          open(OUT,">> $read_out") or die "Error writing to $read_out\n";
          print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tTransplicing";
          close OUT;        
        }
        else
        {
          open(OUT,">> $read_out") or die "Error writing to $read_out\n";
          print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tInvestigate";
          close OUT;       
        } 
      }

      elsif($gene1mel eq "*" && $gene1sec ne "*") # mate1 aligned to sec but not mel
      {
        if($gene2sec ne "*" && $gene1sec eq $gene2sec && $sec2mel{$gene1sec}) #7 & 8
        {
          $geneHash{$sec2mel{$gene1sec}}{$species2} += 1; 
          if($exon1sec == $exon2sec){$exonHash{$sec2mel{$gene1sec}}{$exon1sec}{$species2} += 1;} 
          else{$exonHash{$sec2mel{$gene1sec}}{$exon1sec}{$species2} += 1; $exonHash{$sec2mel{$gene1sec}}{$exon2sec}{$species2} += 1;}  
          open(OUT,">> $read_out") or die "Error writing to $read_out\n";
          print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\t$species2";
          close OUT;        
        } 
        elsif($gene2mel eq "*" && $gene2sec eq "*" && $sec2mel{$gene1sec}) #9
        {
          $geneHash{$sec2mel{$gene1sec}}{$species2} += 1; 
          $exonHash{$sec2mel{$gene1sec}}{$exon1sec}{$species2} += 1; 
          open(OUT,">> $read_out") or die "Error writing to $read_out\n";
          print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\t$species2";
          close OUT; 
        }
        elsif($gene2mel ne "*" && $gene2sec eq "*") #15
        {
          open(OUT,">> $read_out") or die "Error writing to $read_out\n";
          print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tTransplicing";
          close OUT;       
        } 
        else
        {
          open(OUT,">> $read_out") or die "Error writing to $read_out\n";
          print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tInvestigate";
          close OUT;       
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
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\t$species1";
            close OUT; 
          }     
          elsif($mel2CG{$gene2mel})
	  {
            $geneHash{$mel2CG{$gene2mel}}{$species1} += 1; 
            $exonHash{$mel2CG{$gene2mel}}{$exon2mel}{$species1} += 1;   
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\t$species1";
            close OUT; 
          }
          else
          {
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tError5";
            close OUT;       
          }
        } 
        elsif($gene2mel eq "*" && $gene2sec ne "*" && $sec2mel{$gene2sec}) #10
        {
          $geneHash{$sec2mel{$gene2sec}}{$species2} += 1; 
          $exonHash{$sec2mel{$gene2sec}}{$exon2sec}{$species2} += 1;
          open(OUT,">> $read_out") or die "Error writing to $read_out\n";
          print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\t$species2";
          close OUT; 
        }
        elsif($gene2mel ne "*" && $gene2sec ne "*" && $sec2mel{$gene2sec}) #12
        {
          if($gene2mel eq $sec2mel{$gene2sec})
	  {
            $geneHash{$gene2mel}{b} += 1; 
            $exonHash{$gene2mel}{$exon2mel}{b} += 1;   
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tBoth";
            close OUT;    
	  }
          elsif($mel2CG{$gene2mel} eq $sec2mel{$gene2sec})
	  {
            $geneHash{$mel2CG{$gene2mel}}{b} += 1; 
            $exonHash{$mel2CG{$gene2mel}}{$exon2mel}{b} += 1;   
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tBoth";
            close OUT;    
	  }
          else
          {
            open(OUT,">> $read_out") or die "Error writing to $read_out\n";
            print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tError10";
            close OUT;       
          }
        } 
        elsif($gene2mel eq "*" && $gene2sec eq "*") #16
        {
          open(OUT,">> $read_out") or die "Error writing to $read_out\n";
          print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tNA";
          close OUT;       
        }
        else
        {
          open(OUT,">> $read_out") or die "Error writing to $read_out\n";
          print OUT "\n$read1mel\t$seq1mel\t$seq2mel\t$align1mel\t$align1sec\t$align2mel\t$align2sec\tInvestigate";
          close OUT;       
        }
      }
    }       
  }
}
close SAM1;
close SAM2;
close SAM3;
close SAM4;

my$exon_out = $ARGV[8];
open(OUT,">> $exon_out") or die "Error writing to $exon_out\n";
print OUT "gene\texon\t$species1\t$species2\tBoth\tAdj$species1\tAdj$species2";
close OUT;

my$madj = 0;
my$sadj = 0;

foreach my$k1(sort keys %exonHash)
{
  foreach my$k2 (sort keys %{$exonHash{$k1}})
  {
    if($exonHash{$k1}{$k2}{$species1} && $exonHash{$k1}{$k2}{$species2} && $exonHash{$k1}{$k2}{b})
    {
      $madj = $exonHash{$k1}{$k2}{$species1} + ($exonHash{$k1}{$k2}{$species1} / ($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2})) * $exonHash{$k1}{$k2}{b};  
      $sadj = $exonHash{$k1}{$k2}{$species2} + ($exonHash{$k1}{$k2}{$species2} / ($exonHash{$k1}{$k2}{$species1} + $exonHash{$k1}{$k2}{$species2})) * $exonHash{$k1}{$k2}{b};

      open(OUT,">> $exon_out") or die "Error writing to $exon_out\n";
      print OUT "\n$k1\t$k2\t$exonHash{$k1}{$k2}{$species1}\t$exonHash{$k1}{$k2}{$species2}\t$exonHash{$k1}{$k2}{b}\t$madj\t$sadj";
      close OUT;
    }
    elsif($exonHash{$k1}{$k2}{$species1} && $exonHash{$k1}{$k2}{$species2})
    {
      open(OUT,">> $exon_out") or die "Error writing to $exon_out\n";
      print OUT "\n$k1\t$k2\t$exonHash{$k1}{$k2}{$species1}\t$exonHash{$k1}{$k2}{$species2}\t0\t$exonHash{$k1}{$k2}{$species1}\t$exonHash{$k1}{$k2}{$species2}";
      close OUT;
    }
    elsif($exonHash{$k1}{$k2}{$species1} && $exonHash{$k1}{$k2}{b})
    {
      $madj = $exonHash{$k1}{$k2}{$species1} + 0.5 * $exonHash{$k1}{$k2}{b};  

      open(OUT,">> $exon_out") or die "Error writing to $exon_out\n";
      print OUT "\n$k1\t$k2\t$exonHash{$k1}{$k2}{$species1}\t0\t$exonHash{$k1}{$k2}{b}\t$madj\t0";
      close OUT;
    }
    elsif($exonHash{$k1}{$k2}{$species2} && $exonHash{$k1}{$k2}{b})
    {
      $sadj = $exonHash{$k1}{$k2}{$species2} + 0.5 * $exonHash{$k1}{$k2}{b};

      open(OUT,">> $exon_out") or die "Error writing to $exon_out\n";
      print OUT "\n$k1\t$k2\t0\t$exonHash{$k1}{$k2}{$species2}\t$exonHash{$k1}{$k2}{b}\t0\t$sadj";
      close OUT;
    }
    elsif($exonHash{$k1}{$k2}{$species1})
    {
      open(OUT,">> $exon_out") or die "Error writing to $exon_out\n";
      print OUT "\n$k1\t$k2\t$exonHash{$k1}{$k2}{$species1}\t0\t0\t$exonHash{$k1}{$k2}{$species1}\t0";
      close OUT;
    }
    elsif($exonHash{$k1}{$k2}{$species2})
    {
      open(OUT,">> $exon_out") or die "Error writing to $exon_out\n";
      print OUT "\n$k1\t$k2\t0\t$exonHash{$k1}{$k2}{$species2}\t0\t0\t$exonHash{$k1}{$k2}{$species2}";
      close OUT;
    }
    elsif($exonHash{$k1}{$k2}{b})
    {
      open(OUT,">> $exon_out") or die "Error writing to $exon_out\n";
      print OUT "\n$k1\t$k2\t0\t0\t$exonHash{$k1}{$k2}{b}\t0\t0";
      close OUT;
    }
    else
    {
      open(OUT,">> $exon_out") or die "Error writing to $exon_out\n";
      print OUT "\n$k1\t$k2\t0\t0\t0";
      close OUT;
    }
    $madj = 0;
    $sadj = 0;
  }
}

my$gene_out = $ARGV[9];
open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
print OUT "gene\t$species1\t$species2\tBoth\tAdj$species1\tAdj$species2";
close OUT;

$madj = 0;
$sadj = 0;

foreach(sort keys %geneHash)
{
  if($geneHash{$_}{$species1} && $geneHash{$_}{$species2} && $geneHash{$_}{b})
  {
    $madj = $geneHash{$_}{$species1} + ($geneHash{$_}{$species1} / ($geneHash{$_}{$species1} + $geneHash{$_}{$species2})) * $geneHash{$_}{b};  
    $sadj = $geneHash{$_}{$species2} + ($geneHash{$_}{$species2} / ($geneHash{$_}{$species1} + $geneHash{$_}{$species2})) * $geneHash{$_}{b};

    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}\t$geneHash{$_}{b}\t$madj\t$sadj";
    close OUT;
  }
  elsif($geneHash{$_}{$species1} && $geneHash{$_}{$species2})
  {
    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$_\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}\t0\t$geneHash{$_}{$species1}\t$geneHash{$_}{$species2}";
    close OUT;
  }
  elsif($geneHash{$_}{$species1} && $geneHash{$_}{b})
  {
    $madj = $geneHash{$_}{$species1} + 0.5 * $geneHash{$_}{b};  
  
    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$_\t$geneHash{$_}{$species1}\t0\t$geneHash{$_}{b}\t$madj\t0";
    close OUT;
  }
  elsif($geneHash{$_}{$species2} && $geneHash{$_}{b})
  {
    $sadj = $geneHash{$_}{$species2} + 0.5 * $geneHash{$_}{b};

    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$_\t0\t$geneHash{$_}{$species2}\t$geneHash{$_}{b}\t0\t$sadj";
    close OUT;
  }
  elsif($geneHash{$_}{$species1})
  {
    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$_\t$geneHash{$_}{$species1}\t0\t0\t$geneHash{$_}{$species1}\t0";
    close OUT;
  }
  elsif($geneHash{$_}{$species2})
  {
    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$_\t0\t$geneHash{$_}{$species2}\t0\t0\t$geneHash{$_}{$species2}";
    close OUT;
  }
  elsif($geneHash{$_}{b})
  {
    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$_\t0\t0\t$geneHash{$_}{b}\t0\t0";
    close OUT;
  }
  else
  {
    open(OUT,">> $gene_out") or die "Error writing to $gene_out\n";
    print OUT "\n$_\t0\t0\t0";
    close OUT;
  }
  $madj = 0;
  $sadj = 0;
}
