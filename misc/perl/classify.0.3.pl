#!/usr/bin/perl

###################################################################################
#
# 06/09/2010
#
# classify.0.3.pl
#
# Purpose: combine 2 .bed files (lifted) to elucidate allele-specific expression 
# 
# Input: two separate alignment files (1 for each species)
#
# Output: a list of genes and their respective match and mismatch counts, along with a 
#         summary of all the reads 
#
# E.g. perl classify.0.2.pl <species1> <species2> <bedmate1species1> <bedmate2species1> <bedmate1species2> <bedmate2species2> <readout> <exonout> <geneout>
#      <species1>         ==> indicate name of 1st species aligned to (i.e. Dmel, Dsec, etc.)
#      <species2>         ==> indicate name of 2nd species aligned to
#      <bedspecies1> ==> file with mate 1 alignments to 1st species
#      <bedspecies2> ==> file with mate 1 alignments to 2nd species
#      <readout>          ==> file to send reads for output
#      <exonout>          ==> file to send exon counts for output
#      <geneout>          ==> file to send gene counts for output
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
my@elements;
my%geneHash;
my%exonHash;
my%mate1spec1;
my%mate2spec1;
my%mate1spec2;
my%mate2spec2;
my%readList;

my$species1 = $ARGV[0]; # Dmel
my$species2 = $ARGV[1]; # Dsec
my$bedspecies1 = $ARGV[2];
my$bedspecies2 = $ARGV[3];
my$read_out = $ARGV[4]; 
my$exon_out = $ARGV[5];
my$gene_out = $ARGV[6];

# bedfile species 1
open(BED,"$bedspecies1") or die "Can't open $bedspecies1 for reading!";
LINE:while(<BED>) 
{
  chomp;
  $line = $_;  
  #print "new line\n";
  if($line =~ /^\@.+/){next LINE;}
  elsif($line =~ /HWI.+\/1/)
  { 
    @elements = split(/\s/,$line);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene = $elements[0];
      $exon = $elements[1] . "," . $elements[2];
      $mate1spec1{$read} = [$gene,$exon]; 
      $readList{$read} = 1;
    }
  }
  elsif($line =~ /HWI.+\/2/)
  {
    @elements = split(/\s/,$line);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene = $elements[0];
      $exon = $elements[1] . "," . $elements[2];
      $mate2spec1{$read} = [$gene,$exon]; 
      $readList{$read} = 1;
    }
  }
} 
close BED;   

print "\nDone making read hashes for 1st species!\n";

# bedfile species 2
open(BED,"$bedspecies2") or die "Can't open $bedspecies2 for reading!";
LINE:while(<BED>) 
{
  chomp;
  $line = $_;  
  #print "new line\n";
  if($line =~ /^\@.+/){next LINE;}
  elsif($line =~ /HWI.+\/1/)
  { 
    @elements = split(/\s/,$line);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene = $elements[0];
      $exon = $elements[1] . "," . $elements[2];
      $mate1spec2{$read} = [$gene,$exon]; 
      $readList{$read} = 1;
    }
  }
  elsif($line =~ /HWI.+\/2/)
  {
    @elements = split(/\s/,$line);
    if($elements[3] =~ /^(HWI.+)\/\d{1}/){$read = $1;}
    if($elements[0] =~ /^C\w{1}\d+/)
    {
      $gene = $elements[0];
      $exon = $elements[1] . "," . $elements[2];
      $mate2spec2{$read} = [$gene,$exon]; 
      $readList{$read} = 1;
    }
  }
}
close BED;

print "\nDone making read hashes for 2nd species!\n";

##### start counting #####
#there should be 16 possibilities for alignments of both mates to both species (4 files, 2 possibilities each, 2^4 = 16)
#these can be taken care of on a case by case basis

open(OUT1,">$read_out") or die "Error writing to $read_out\n";
print OUT1 "read\tgene1$species1\tgene1$species2\tgene2$species1\tgene2$species2\tcall"; 

foreach $read (keys %readList)   
{
  if($mate1spec1{$read} && $mate1spec2{$read}) # mate1 aligned to both genomes
  {
    if($mate2spec1{$read} && $mate2spec2{$read} && $mate1spec1{$read}[0] eq $mate2spec1{$read}[0] && $mate1spec2{$read}[0] eq $mate2spec2{$read}[0]) #11
    {
      if($mate1spec1{$read}[0] eq $mate1spec2{$read}[0])
      {
        $geneHash{$mate1spec1{$read}[0]}{b} += 1; 
        if($mate1spec1{$read}[1] eq $mate2spec1{$read}[1]){$exonHash{$mate1spec1{$read}[0]}{$mate1spec1{$read}[1]}{b} += 1;} 
        else{$exonHash{$mate1spec1{$read}[0]}{$mate1spec1{$read}[1]}{b} += 1; $exonHash{$mate1spec1{$read}[0]}{$mate2spec1{$read}[1]}{b} += 1;}     
        print OUT1 "\n$read\t$mate1spec1{$read}[0]\t$mate1spec2{$read}[0]\t$mate2spec1{$read}[0]\t$mate2spec2{$read}[0]\tBoth";
      } 
      else{print OUT1 "\n$read\t$mate1spec1{$read}[0]\t$mate1spec2{$read}[0]\t$mate2spec1{$read}[0]\t$mate2spec2{$read}[0]\tError11";}
    } 
    elsif($mate2spec1{$read}) #1
    {
      if($mate1spec1{$read}[0] eq $mate2spec1{$read}[0])
      {
        $geneHash{$mate1spec1{$read}[0]}{$species1} += 1; 
        if($mate1spec1{$read}[1] eq $mate2spec1{$read}[1]){$exonHash{$mate1spec1{$read}[0]}{$mate1spec1{$read}[1]}{$species1} += 1;} 
        else{$exonHash{$mate1spec1{$read}[0]}{$mate1spec1{$read}[1]}{$species1} += 1; $exonHash{$mate1spec1{$read}[0]}{$mate2spec1{$read}[1]}{$species1} += 1;} 
        print OUT1 "\n$read\t$mate1spec1{$read}[0]\t$mate1spec2{$read}[0]\t$mate2spec1{$read}[0]\t\*\t$species1";
      }
      else{print OUT1 "\n$read\t$mate1spec1{$read}[0]\t$mate1spec2{$read}[0]\t$mate2spec1{$read}[0]\t\*\tError1";}
    }
    elsif($mate2spec2{$read}) #6
    {
      if($mate1spec2{$read}[0] eq $mate2spec2{$read}[0])
      {
        $geneHash{$mate1spec2{$read}[0]}{$species2} += 1; 
        if($mate1spec2{$read}[1] eq $mate2spec2{$read}[1]){$exonHash{$mate1spec2{$read}[0]}{$mate1spec2{$read}[1]}{$species2} += 1;} 
        else{$exonHash{$mate1spec2{$read}[0]}{$mate1spec2{$read}[1]}{$species2} += 1; $exonHash{$mate1spec2{$read}[0]}{$mate2spec2{$read}[1]}{$species2} += 1;} 
        print OUT1 "\n$read\t$mate1spec1{$read}[0]\t$mate1spec2{$read}[0]\t\*\t$mate2spec2{$read}[0]\t$species2";
      }
      else{print OUT1 "\n$read\t$mate1spec1{$read}[0]\t$mate1spec2{$read}[0]\t\*\t$mate2spec2{$read}[0]\tError6";}
    }
    else #13
    {
      if($mate1spec1{$read}[0] eq $mate1spec2{$read}[0] && $mate1spec1{$read}[1] eq $mate1spec2{$read}[1])
      {
        $geneHash{$mate1spec1{$read}[0]}{b} += 1; 
        $exonHash{$mate1spec1{$read}[0]}{$mate1spec1{$read}[1]}{b} += 1;    
        print OUT1 "\n$read\t$mate1spec1{$read}[0]\t$mate1spec2{$read}[0]\t\*\t\*\tBoth";
      }       
      else{print OUT1 "\n$read\t$mate1spec1{$read}[0]\t$mate1spec2{$read}[0]\t\*\t\*\tError13";}
    } 
  }

  elsif($mate1spec1{$read}) # mate1 aligned to mel but not sec
  {
    if($mate2spec1{$read} && $mate2spec2{$read}) #2
    {
      if($mate1spec1{$read}[0] eq $mate2spec1{$read}[0])
      {
        $geneHash{$mate1spec1{$read}[0]}{$species1} += 1; 
        if($mate1spec1{$read}[1] eq $mate2spec1{$read}[1]){$exonHash{$mate1spec1{$read}[0]}{$mate1spec1{$read}[1]}{$species1} += 1;} 
        else{$exonHash{$mate1spec1{$read}[0]}{$mate1spec1{$read}[1]}{$species1} += 1; $exonHash{$mate1spec1{$read}[0]}{$mate2spec1{$read}[1]}{$species1} += 1;} 
        print OUT1 "\n$read\t$mate1spec1{$read}[0]\t\*\t$mate2spec1{$read}[0]\t$mate2spec2{$read}[0]\t$species1";
      }
      else{print OUT1 "\n$read\t$mate1spec1{$read}[0]\t\*\t$mate2spec1{$read}[0]\t$mate2spec2{$read}[0]\tError2";}
    } 
    elsif($mate2spec1{$read}) #3
    {
      if($mate1spec1{$read}[0] eq $mate2spec1{$read}[0])
      {
        $geneHash{$mate1spec1{$read}[0]}{$species1} += 1; 
        if($mate1spec1{$read}[1] eq $mate2spec1{$read}[1]){$exonHash{$mate1spec1{$read}[0]}{$mate1spec1{$read}[1]}{$species1} += 1;} 
        else{$exonHash{$mate1spec1{$read}[0]}{$mate1spec1{$read}[1]}{$species1} += 1; $exonHash{$mate1spec1{$read}[0]}{$mate2spec1{$read}[1]}{$species1} += 1;} 
        print OUT1 "\n$read\t$mate1spec1{$read}[0]\t\*\t$mate2spec1{$read}[0]\t\*\t$species1";
      }
      else{print OUT1 "\n$read\t$mate1spec1{$read}[0]\t\*\t$mate2spec1{$read}[0]\t\*\tError3";}
    }
    elsif($mate2spec2{$read}) #14
    {
      print OUT1 "\n$read\t$mate1spec1{$read}[0]\t\*\t\*\t$mate2spec2{$read}[0]\tTransplicing";
    }
    else #4
    {
      $geneHash{$mate1spec1{$read}[0]}{$species1} += 1; 
      $exonHash{$mate1spec1{$read}[0]}{$mate1spec1{$read}[1]}{$species1} += 1; 
      print OUT1 "\n$read\t$mate1spec1{$read}[0]\t\*\t\*\t\*\t$species1";
    }       
  }

  elsif($mate1spec2{$read}) # mate1 aligned to sec but not mel
  {
    if($mate2spec1{$read} && $mate2spec2{$read}) #7
    {
      if($mate1spec2{$read}[0] eq $mate2spec2{$read}[0])
      {
        $geneHash{$mate1spec2{$read}[0]}{$species2} += 1; 
        if($mate1spec2{$read}[1] eq $mate2spec2{$read}[1]){$exonHash{$mate1spec2{$read}[0]}{$mate1spec2{$read}[1]}{$species2} += 1;} 
        else{$exonHash{$mate1spec2{$read}[0]}{$mate1spec2{$read}[1]}{$species2} += 1; $exonHash{$mate1spec2{$read}[0]}{$mate2spec2{$read}[1]}{$species2} += 1;}  
        print OUT1 "\n$read\t\*\t$mate1spec2{$read}[0]\t$mate2spec1{$read}[0]\t$mate2spec2{$read}[0]\t$species2";
      }
      else{print OUT1 "\n$read\t\*\t$mate1spec2{$read}[0]\t$mate2spec1{$read}[0]\t$mate2spec2{$read}[0]\tError7";}
    } 
    elsif($mate2spec2{$read}) #8
    {
      if($mate1spec2{$read}[0] eq $mate2spec2{$read}[0])
      {
        $geneHash{$mate1spec2{$read}[0]}{$species2} += 1; 
        if($mate1spec2{$read}[1] eq $mate2spec2{$read}[1]){$exonHash{$mate1spec2{$read}[0]}{$mate1spec2{$read}[1]}{$species2} += 1;} 
        else{$exonHash{$mate1spec2{$read}[0]}{$mate1spec2{$read}[1]}{$species2} += 1; $exonHash{$mate1spec2{$read}[0]}{$mate2spec2{$read}[1]}{$species2} += 1;}  
        print OUT1 "\n$read\t\*\t$mate1spec2{$read}[0]\t\*\t$mate2spec2{$read}[0]\t$species2";
      }
      else{print OUT1 "\n$read\t\*\t$mate1spec2{$read}[0]\t\*\t$mate2spec2{$read}[0]\tError8";}
    }
    elsif($mate2spec1{$read}) #15
    {
      print OUT1 "\n$read\t\*\t$mate1spec2{$read}[0]\t$mate2spec1{$read}[0]\t\*\tTransplicing";
    } 
    else #9
    {
      $geneHash{$mate1spec2{$read}[0]}{$species2} += 1; 
      $exonHash{$mate1spec2{$read}[0]}{$mate1spec2{$read}[1]}{$species2} += 1; 
      print OUT1 "\n$read\t\*\t$mate1spec2{$read}[0]\t\*\t\*\t$species2";        
    }
  }

  else # mate1 aligned to nothing
  {
    if($mate2spec1{$read} && $mate2spec2{$read}) #12
    {
      if($mate2spec1{$read}[0] eq $mate2spec2{$read}[0] && $mate2spec1{$read}[1] eq $mate2spec2{$read}[1])
      {
        $geneHash{$mate2spec1{$read}[0]}{b} += 1; 
        $exonHash{$mate2spec1{$read}[0]}{$mate2spec1{$read}[1]}{b} += 1;   
        print OUT1 "\n$read\t\*\t\*\t$mate2spec1{$read}[0]\t$mate2spec2{$read}[0]\tBoth";
      }
      else{print OUT1 "\n$read\t\*\t\*\t$mate2spec1{$read}[0]\t$mate2spec2{$read}[0]\tError12";}
    }
    elsif($mate2spec1{$read}) #5
    {
      $geneHash{$mate2spec1{$read}[0]}{$species1} += 1; 
      $exonHash{$mate2spec1{$read}[0]}{$mate2spec1{$read}[1]}{$species1} += 1;   
      print OUT1 "\n$read\t\*\t\*\t$mate2spec1{$read}[0]\t\*\t$species1";              
    } 
    elsif($mate2spec2{$read}) #10
    {
      $geneHash{$mate2spec2{$read}[0]}{$species2} += 1; 
      $exonHash{$mate2spec2{$read}[0]}{$mate2spec2{$read}[1]}{$species2} += 1;
      print OUT1 "\n$read\t\*\t\*\t\*\t$mate2spec2{$read}[0]\t$species2";
    }
    else #16
    {
      print OUT1 "\n$read\t\*\t\*\t\*\t\*\tNA";
    }       
  }
}
close OUT1;

open(OUT2,"> $exon_out") or die "Error writing to $exon_out\n";
print OUT2 "gene\texon\t$species1\t$species2\tBoth\tAdj$species1\tAdj$species2";

my$madj = 0;
my$sadj = 0;

foreach my$k1 (keys %exonHash)
{
  foreach my$k2 (keys %{$exonHash{$k1}})
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
      $madj = $exonHash{$k1}{$k2}{$species1} + 0.5*$exonHash{$k1}{$k2}{b};  
      print OUT2 "\n$k1\t$k2\t$exonHash{$k1}{$k2}{$species1}\t0\t$exonHash{$k1}{$k2}{b}\t$madj\t0";
    }
    elsif($exonHash{$k1}{$k2}{$species2} && $exonHash{$k1}{$k2}{b})
    {
      $sadj = $exonHash{$k1}{$k2}{$species2} + 0.5*$exonHash{$k1}{$k2}{b};
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

open(OUT3,"> $gene_out") or die "Error writing to $gene_out\n";
print OUT3 "gene\t$species1\t$species2\tBoth\tAdj$species1\tAdj$species2";

$madj = 0;
$sadj = 0;

foreach(keys %geneHash)
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
    $madj = $geneHash{$_}{$species1} + 0.5*$geneHash{$_}{b};  
    print OUT3 "\n$_\t$geneHash{$_}{$species1}\t0\t$geneHash{$_}{b}\t$madj\t0";
  }
  elsif($geneHash{$_}{$species2} && $geneHash{$_}{b})
  {
    $sadj = $geneHash{$_}{$species2} + 0.5*$geneHash{$_}{b};
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
}
