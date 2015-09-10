#!/usr/bin/perl
#
# exon_hitter_read.pl (originally exon_iso_gene_hits.pl) 
# duff 10/30/2008
#
# Read in Joel's filtered/post-processed bowtie alignments
# and jodi's ucsc/flybase constitutive exons list.
# Find reads that map to exons and write read,exon, and gene info to file

# NOTE: GENES AND EXONS ARE DEFINED BY CONSTITUTIVE EXONS ONLY!!!!
# ALSO, THIS TIME I'M USING JODI'S NEW CONSTITUTIVE EXON LIST!!!

use strict;
use FileHandle;



my @READS_FILENAME;
$READS_FILENAME[1] = "Hyb_mel_nogap_unq.txt";
$READS_FILENAME[2] = "Hyb_sec_lift_nogap_unq.txt";
$READS_FILENAME[3] = "Mix_mel_nogap_unq.txt";
$READS_FILENAME[4] = "Mix_sec_lift_nogap_unq.txt";



open(CONS_EXONS, "<constitutive_regions_final.bed") or die("can't open all_cons_exons_inc_chrM_fix.txt"); #UCSC/FLYBASE genes table
<CONS_EXONS>; 

my @f;
my @g;
my @exon_start;
my @exon_end;

my $line;


my %exon_to_gene;  #{chromo}{start}{end}->[gene]    *****Assumes that exons (defined by <chrmo,start,end> belong to only one gene.
                                                     # ***** CORRECTED TO ALLOW FOR LIST OF GENES FOR EACH EXON
my %gene_cross_section; #{gene}->cross-section
my $exon_cross_section;

my $chromo;
my $start;
my $end;
my $strand;

my $gene;


my %unique_exons_start;
my %unique_exons_end;

my $estart;
my $eend;

my $rpkm;

#-------------

my $read_id;

my %gene_list;



#PROCESS EXON & GENE INFORMATION  (do this once)

printf("Processing exon and gene coordinate information...");

while ($line = <CONS_EXONS>){
  chomp($line);
  @g = split /\s+/,$line;  
  $gene = $g[3];
  $chromo = $g[0];
  $gene_list{$gene}=1;
  $estart = $g[1]+1;
  $eend   = $g[2];
  push(@{$exon_to_gene{$chromo}{$estart}{$eend}}, $gene);            #*****Assumes that exons (defined by <chrmo,start,end> belong to only one gene.

}


#unique exons to equivalent arrays and compute gene cross-section (sum of constitutive exon lengths):

foreach my $chromo (keys %exon_to_gene){                  # for each chromo
                                                             # sort exons by starts
  foreach my $estart (keys %{$exon_to_gene{$chromo}}){
    foreach my $eend (keys %{$exon_to_gene{$chromo}{$estart}}){
      push(@{$unique_exons_start{$chromo}},$estart);
      push(@{$unique_exons_end{$chromo}},$eend); 
      for(my $i=0;$i<=$#{@{$exon_to_gene{$chromo}{$estart}{$eend}}};$i++){    #loop over all genes that have this exon
        $gene_cross_section{$exon_to_gene{$chromo}{$estart}{$eend}[$i]} += $eend - $estart +1;   #gene cross-section here is determined by cons exons ONLY!!!!
      }
    }
  }
                                        #sort by exon start and mirror the order in the end array
  @{$unique_exons_start{$chromo}} = @{$unique_exons_start{$chromo}}[
                              my @idx = sort {$unique_exons_start{$chromo}[$a] <=> $unique_exons_start{$chromo}[$b]} 0..$#{@{$unique_exons_start{$chromo}}}];

  @{$unique_exons_end{$chromo}} = @{$unique_exons_end{$chromo}}[@idx];
  


printf("Finished with processing exon/gene info....\n");

}

foreach $gene (keys %gene_list){
  if ($gene_cross_section{$gene} == 0){
    printf("$gene has cross section 0 \n");
  }
}

  
#-----------------------------------------   READS  ----------------------------------------------


for(my $ifile=1; $ifile<=4; $ifile++){          #LOOP OVER READ FILES


my %read_read_id;
my %read_start;
my %read_end;
my %read_strand;

my %gene_hit_count;            #these are defined by con exons
my %exon_hit_count;

my $total_exon_aligned_reads = 0;

open(READS, "<$READS_FILENAME[$ifile]") or die("can't open $READS_FILENAME[$ifile]");



my $OUTFILE = $READS_FILENAME[$ifile];
$OUTFILE =~ s/.txt//;
my $CONS_EXON_READS = $OUTFILE."_cons_exon_reads_Jan_new.txt";
my $CONS_EXON_HIT_COUNTS = $OUTFILE."_cons_exon_hit_counts_Jan_new.txt";
my $CONS_GENE_HIT_COUNTS = $OUTFILE."_cons_gene_hit_counts_Jan_new.txt";

my $CONS_EXON_READS = new FileHandle ">$CONS_EXON_READS" or die "can't open $CONS_EXON_READS";   #reads that hit cons exons

my $CONS_EXON_HITS = new FileHandle ">$CONS_EXON_HIT_COUNTS" or die "can't open $CONS_EXON_HIT_COUNTS";   #output file of exon hit counts

my $CONS_GENE_HITS = new FileHandle ">$CONS_GENE_HIT_COUNTS" or die "can't open $CONS_GENE_HIT_COUNTS";   #output file of gene hit counts


printf("-----------------------reading reads from $READS_FILENAME[$ifile]...\n");
my $linecount=0;
while($line= <READS>){
  chomp($line);
  @f = split /\s+/,$line;
  $chromo  = $f[2];
  $start   = $f[3] + 1;
  $end     = $f[4];
  $read_id = $f[5];
  $strand  = $f[7];

  push(@{$read_read_id{$chromo}}, $read_id);
  push(@{$read_start{$chromo}}, $start);
  push(@{$read_end{$chromo}}, $end);
  push(@{$read_strand{$chromo}}, $strand);

  $linecount++;
  if ($linecount%1000 == 0){
    printf("$linecount read: $chromo $start $end $strand \n");
  }
}

  

#reads sorted by start for each chromo

printf("sorting reads... \n");

foreach $chromo (keys %read_start){


  @{$read_start{$chromo}} = @{$read_start{$chromo}}[
                              my @idx = sort {$read_start{$chromo}[$a] <=> $read_start{$chromo}[$b]} 0..$#{@{$read_start{$chromo}}}];

  @{$read_read_id{$chromo}} = @{$read_read_id{$chromo}}[@idx];
  @{$read_end{$chromo}} = @{$read_end{$chromo}}[@idx];
  @{$read_strand{$chromo}} = @{$read_strand{$chromo}}[@idx];

}
 


#FOR EACH READ, FIND ITS EXON COORDINATES AND GENE NAME IF IT HITS

printf("processing reads for exon hits ....\n");

$linecount=0;


#---------



foreach $chromo (keys %read_start){

  my $order_index = 0;

  for (my $i=0; $i <= $#{@{$read_start{$chromo}}}; $i++){

     $start = $read_start{$chromo}[$i];
     $end   = $read_end{$chromo}[$i];
     $read_id = $read_read_id{$chromo}[$i];
#printf(" READ $start $end \n");
     $estart = $unique_exons_start{$chromo}[$order_index];   #initial value
     my $hit_flag=0;
     while ($estart <= $start && $order_index <= $#{@{$unique_exons_start{$chromo}}} && ($hit_flag != 1)){
       $estart = $unique_exons_start{$chromo}[$order_index];
       $eend = $unique_exons_end{$chromo}[$order_index];
#printf("  EXON $estart $eend \n");
       if ($start >= $estart && $end <= $eend){
         $CONS_EXON_READS->printf("$read_id\t$chromo\t$start\t$end\t$estart\t$eend\t"); 
         my $gene_string;         
         for(my $i=0;$i<=$#{@{$exon_to_gene{$chromo}{$estart}{$eend}}};$i++){    #loop over all genes that have this exon
           $gene = $exon_to_gene{$chromo}{$estart}{$eend}[$i];
           $gene_hit_count{$gene}++;     
           $gene_string = $gene_string."_".$gene;       
	 }
         $gene_string = substr($gene_string,1,length($gene_string)-1);
         $CONS_EXON_READS->printf("$gene_string\n");
         $exon_hit_count{$chromo}{$estart}{$eend}++;
	 $hit_flag = 1;
         $total_exon_aligned_reads++;
       }
       if ($hit_flag != 1){
         $order_index++;
         $estart = $unique_exons_start{$chromo}[$order_index];      #This line was added 11/6/08        
       }
     }  #next exon


   } #next read

} #next chromo



close(READS);

printf("total_exon_aligned_reads  $total_exon_aligned_reads \n");

printf("DONE.\n");


foreach $gene (keys %gene_list){
    if (!exists $gene_hit_count{$gene}){
	$gene_hit_count{$gene} = 0;
    }
    $rpkm = $gene_hit_count{$gene}/($gene_cross_section{$gene}/1000)/($total_exon_aligned_reads/1000000);
    $CONS_GENE_HITS->printf("$gene\t$gene_hit_count{$gene}\t$gene_cross_section{$gene}\t$rpkm\n");
    
}


foreach my $chromo (keys %exon_to_gene){                  # for each exon
   foreach my $estart (keys %{$exon_to_gene{$chromo}}){
    foreach my $eend (keys %{$exon_to_gene{$chromo}{$estart}}){

      if (!exists $exon_hit_count{$chromo}{$estart}{$eend} ){
	$exon_hit_count{$chromo}{$estart}{$eend} = 0;
      }
      $exon_cross_section = $eend - $estart + 1;
      $rpkm = $exon_hit_count{$chromo}{$estart}{$eend}/($exon_cross_section/1000)/($total_exon_aligned_reads/1000000);   


      my $gene_string;         
      for(my $i=0;$i<=$#{@{$exon_to_gene{$chromo}{$estart}{$eend}}};$i++){    #loop over all genes that have this exon
         $gene = $exon_to_gene{$chromo}{$estart}{$eend}[$i];
         $gene_string = $gene_string."_".$gene;       
      }
      $gene_string = substr($gene_string,1,length($gene_string)-1);



   
      $CONS_EXON_HITS->printf("$gene_string\t$chromo\t$estart\t$eend\t$exon_hit_count{$chromo}{$estart}{$eend}\t$exon_cross_section\t$rpkm\n");



    }
  }
 }



$CONS_EXON_READS->close;

$CONS_EXON_HITS->close; 

$CONS_GENE_HITS->close; 




} #next file














