#!/usr/bin/perl

#######################
# 09/15/2009
#
# BLAST2contigs.pl
#
# Purpose: using short sequence reads, build a list of top genes for each
#          sequence read using a BLAST (basic local alignment search tool) query
# 
# Input: each short sequence will be BLASTed compared to the following species:
#        1). D. melanogaster
#        2). D. sechellia
#	 3). D. simulans
#	 4). D. yakuba
#
# Output: collect best hit information for each species including chromosome,
#         local alignment, and gene
#######################

use strict;
use warnings;

my$start = time;

my$data_dir = $ARGV[0];

# user input in case the references change
my$dmelref = $ARGV[1];
my$dsecref = $ARGV[2];
my$dsimref = $ARGV[3];
my$dyakref = $ARGV[4];

# build query databases
#system "makeblastdb -in \$dmelref -dbtype nucl";
#system "makeblastdb -in \$dsecref -dbtype nucl";
#system "makeblastdb -in \$dsimref -dbtype nucl";
#system "makeblastdb -in \$dyakref -dbtype nucl";

opendir(DIR,$data_dir) || die "\nCan't open directory $data_dir for reading\n";
#my@files = grep{/^s_/} grep{/.fa$/} readdir(DIR);
my@files = grep{/^trial1/} grep{/.fa$/} readdir(DIR);
close(DIR);

#foreach(@files){print "\n$_\n\n";}

open (OUT,">> $data_dir/greatestHits.txt") or die "Error writing to $data_dir/greatestHits.txt\n";
print OUT "info\tseq\tdmelBHgene\tdsecBHgene\tdsimBHgene\tdyakBHgene";

foreach(@files)
{   
  system "./blastall -p blastn -d $dmelref -i $_ -m 8 -b 1 > dmelBLAST.txt";
  system "./blastall -p blastn -d $dsecref -i $_ -m 8 -b 1 > dsecBLAST.txt";
  system "./blastall -p blastn -d $dsimref -i $_ -m 8 -b 1 > dsimBLAST.txt";
  system "./blastall -p blastn -d $dyakref -i $_ -m 8 -b 1 > dyakBLAST.txt";
                          
  open(IN,"$data_dir\/$_") or die "\nError opening $data_dir\/$_\n";
  my$i;
  my%seqs;
  my%infomel;
  my%infosec;
  my%infosim;
  my%infoyak;
  my$key;  

  while(<IN>)
  {
    chomp $_;  
    if(/>(.+)\s+/){$key = $1;}
    else{$seqs{$key}=$_;}
  }
  close IN;

  #foreach $key (keys %seqs) {print "\nKey: $key Seq: $seqs{$key}"};

  open(MEL,"$data_dir\/dmelBLAST.txt") or die "\nError opening $data_dir\/dmelBLAST.txt\n";
  while(<MEL>)
  {   
    chomp $_;
    my@array = split("\t",$_);
    $key = $array[0];
    $infomel{$key} = $array[2];
    #print "The key is: {$key} which points to: $infomel{$key}\n";
    #print "It should point to $array[2]\n";
  }
  close MEL;

#   print "seqs\n";
#   foreach $i (keys %seqs)
#   {
#     print "\{$i\}\n";
#   }
#   print "infomel\n";
#   foreach $i (keys %infomel)
#   {
#     print "\{$i\}\n";
#   }

  open(SEC,"$data_dir\/dsecBLAST.txt") or die "\nError opening $data_dir\/dsecBLAST.txt\n";
  while(<SEC>)
  {   
    chomp $_;
    my@array = split("\t",$_);
    $key = $array[0];
    $infosec{$key} = (split(/\\/,$array[2]))[1];  
    #print "The key is: {$key} which points to: $infosec{$key}\n";
    #printf ("It should point to %s\n",(split(/\\/,$array[2]))[1]);
  }
  close SEC;

#   print "seqs\n";
#   foreach $i (keys %seqs)
#   {
#     print "\{$i\}\n";
#   }
#   print "infosec\n";
#   foreach $i (keys %infosec)
#   {
#     print "\{$i\}\n";
#   }

  open(SIM,"$data_dir\/dsimBLAST.txt") or die "\nError opening $data_dir\/dsimBLAST.txt\n";
  while(<SIM>)
  {   
    chomp $_;
    my@array = split("\t",$_);
    $key = $array[0];
    $infosim{$key} = (split(/\\/,$array[2]))[1];  
  }
  close SIM;

# print "seqs\n";
# foreach $i (keys %seqs)
# {
#   print "\{$i\}\n";
# }
# print "infosim\n";
# foreach $i (keys %infosim)
# {
#   print "\{$i\}\n";
# }

  open(YAK,"$data_dir\/dyakBLAST.txt") or die "\nError opening $data_dir\/dyakBLAST.txt\n";
  while(<YAK>)
  {   
    chomp $_;
    my@array = split("\t",$_);
    $key = $array[0];
    $infoyak{$key} = (split(/\\/,$array[2]))[1];  
  }
  close YAK;

 #   print "seqs\n";
 #   foreach $i (keys %seqs)
 #   {
 #     print "\{$i\}\n";
 #   }
 #   print "infoyak\n";
 #   foreach $i (keys %infoyak)
 #   {
 #     print "\{$i\}\n";
 #   }

  foreach $key (keys %infomel)
  {
    print OUT "\n$key\t$seqs{$key}\t$infomel{$key}\t$infosec{$key}\t$infosim{$key}\t$infoyak{$key}";
  }
}
close OUT;
system "rm dmelBLAST.txt dsecBLAST.txt dsimBLAST.txt dyakBLAST.txt";     
printf ("\nTime elapsed: %d\n\n",time-$start);

