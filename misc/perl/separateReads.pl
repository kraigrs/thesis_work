#!/usr/bin/perl

###################################################################################
#
# 08/11/2010
#
# separateReads.pl
#
# Purpose: separate reads that have been classified into respective .fasta[q] files 
# 
# Input: the classified reads from classify.pl, probably called *readout* 
#
# Output: individual fasta[q] files separated into respective files (i.e. separate zhr reads into those correctly 
#         aligned to zhr genome, those aligning to z30 genome, and those aligning equally well to both)
#
# E.g. perl separate.pl <species1> <species2> <readout> <mate1> <mate2>  
#
#      <species1> ==> indicate name of correct species (i.e. Dmel, Dsec, etc.)
#      <species2> ==> indicate name of other species aligned to
#      <readout>  ==> readout file from classify.pl
#      <mate1>    ==> fasta[q] file for mate 1 
#      <mate2>    ==> fasta[q] file for mate 2
#      <output>   ==> prefix for output files
#
###################################################################################

use strict;
use warnings;

main();
sub main
{
  my$species1 = $ARGV[0]; # zhr
  my$species2 = $ARGV[1]; # z30
  my$readout = $ARGV[2];
  my$mate1 = $ARGV[3];
  my$mate2 = $ARGV[4];
  my$output = $ARGV[5]; #prefix for output files

  my$read;
  my$class;
  my@elements;
  my%true;
  my%false;
  my%both;
   
  # open readout file
  open(READS,"$readout") or die "Can't open $readout for reading!";
  LINE:while(<READS>) 
  {
    chomp;  
    if($_ =~ /^read/){next LINE;}
    elsif($_ =~ /^HWI.+/)
    { 
      #print "$_\n";
      @elements = split(/\s/,$_);
      $read = $elements[0];
      $class = $elements[5];
      #print "$read\n";
      if($class eq $species1){$true{$read} = 1;}  
      elsif($class eq $species2){$false{$read} = 1;}
      elsif($class eq "Both"){$both{$read} = 1;}
    }
  }
  close READS;
  print "\nDone making reads hashes!\n";

  open(TRUE1,">$output\.mate1.true.fastq") or die "Error writing to $output\.mate1.true.fastq\n";
  open(FALSE1,">$output\.mate1.false.fastq") or die "Error writing to $output\.mate1.false.fastq\n";
  open(BOTH1,">$output\.mate1.both.fastq") or die "Error writing to $output\.mate1.both.fastq\n";

  open(TRUE2,">$output\.mate2.true.fastq") or die "Error writing to $output\.mate2.true.fastq\n";
  open(FALSE2,">$output\.mate2.false.fastq") or die "Error writing to $output\.mate2.false.fastq\n";
  open(BOTH2,">$output\.mate2.both.fastq") or die "Error writing to $output\.mate2.both.fastq\n";

  my$line1; my$line2; my$read1; my$read2;
  my$foundTrue = 0; my$foundFalse = 0; my$foundBoth = 0;

  # open raw files
  open(MATE1,"$mate1") or die "Can't open $mate1 for reading!";
  open(MATE2,"$mate2") or die "Can't open $mate2 for reading!";
  LINE:while(<MATE1>) 
  {
    chomp;
    $line1 = $_;
    $line2 = <MATE2>;
    chomp $line2;

    if($line1 =~ /^[@+]{1}(HWI.+)\/\d{1}$/){$read1 = $1;} 
    if($line2 =~ /^[@+]{1}(HWI.+)\/\d{1}$/){$read2 = $1;}
     
    if($read1 eq $read2 && $foundTrue == 0 && $foundFalse == 0 && $foundBoth == 0)
    {
      #print "Passed!\n";
      if($true{$read1})
      {
        $foundTrue = 1;
        print TRUE1 "$line1\n";
        print TRUE2 "$line2\n";
        #print "\nThis read is correct\n";
        #print "$line1\n$line2\n";
      }
      elsif($false{$read1})
      {
        $foundFalse = 1;
        print FALSE1 "$line1\n";
        print FALSE2 "$line2\n";
        #print "\nThis read is incorrect\n";
        #print "$line1\n$line2\n";
      }
      elsif($both{$read1})
      {
        $foundBoth = 1;
        print BOTH1 "$line1\n";
        print BOTH2 "$line2\n";
        #print "\nThis read aligns to both\n";
        #print "$line1\n$line2\n";
      }
    }  
    elsif($foundTrue == 1)
    {
      print TRUE1 "$line1\n";
      print TRUE2 "$line2\n";
      $foundTrue = 0;
    }
    elsif($foundFalse == 1)
    {
      print FALSE1 "$line1\n";
      print FALSE2 "$line2\n";
      $foundFalse = 0;
    }
    elsif($foundBoth == 1)
    {
      print BOTH1 "$line1\n";
      print BOTH2 "$line2\n";
      $foundBoth = 0;
    }
  }
  close MATE1; close MATE2; 
  close TRUE1; close TRUE2;
  close FALSE1; close FALSE2;
  close BOTH1; close BOTH2; 
} 
 
