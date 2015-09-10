#!/usr/bin/perl

#######################
# 09/21/2010
#
# convert2SAM.pl
#
# Purpose: converts a long list of .dat files into .sam format and concatenates them
# 
# Input: a directory with a list of .dat files to be converted into .sam (!!! All files must be from same samples !!!)
#
# Output: .sam files corresponding to mates 1 and 2 consisting of all the alignment info originally contained in the .dat files
#######################

use strict;
use warnings;

main ();
sub main 
{ 
  my$data_dir = $ARGV[0];  
  my$prefix = $ARGV[1];
  my$ref1 = $ARGV[2];
  my$ref2 = $ARGV[3];
  
  my@files;

  opendir(DIR,$data_dir) || die "\nCan't open directory $data_dir for reading\n";
  @files = grep{/\.dat\w*$/} readdir(DIR);
  close(DIR);
  #foreach(@files){print "\n$_";}

  foreach(@files){system "MosaikText -in $data_dir/$_ -bed $data_dir/$_\.bed";}
 
  #system "gunzip $data_dir/*.bed.gz";

#  system "cat $data_dir/*_r1_$prefix*$ref1*\.bed > $data_dir/$prefix\.mate1\.$ref1\.mosaik\.bed";
#  system "cat $data_dir/*_r1_$prefix*$ref2*\.bed > $data_dir/$prefix\.mate1\.$ref2\.mosaik\.bed";
#  system "cat $data_dir/*_r2_$prefix*$ref1*\.bed > $data_dir/$prefix\.mate2\.$ref1\.mosaik\.bed";
#  system "cat $data_dir/*_r2_$prefix*$ref2*\.bed > $data_dir/$prefix\.mate2\.$ref2\.mosaik\.bed";

  system "cat $data_dir/$prefix\_r1_$ref1*\.bed > $data_dir/$prefix\.mate1\.$ref1\.mosaik\.bed";
  system "cat $data_dir/$prefix\_r1_$ref2*\.bed > $data_dir/$prefix\.mate1\.$ref2\.mosaik\.bed";
  system "cat $data_dir/$prefix\_r2_$ref1*\.bed > $data_dir/$prefix\.mate2\.$ref1\.mosaik\.bed";
  system "cat $data_dir/$prefix\_r2_$ref2*\.bed > $data_dir/$prefix\.mate2\.$ref2\.mosaik\.bed";
}

