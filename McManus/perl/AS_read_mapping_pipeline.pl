#!/usr/bin/perl
##  Joel McManus 12/30/2008
##  This script runs the whole data analysis pipeline for Hybrid PE sequencing, up to the point of filtering gaps from the results
##  To start, make sure you have the illumina output in fasta format (fastq to fasta conversion scripts are readily available online), 
##  in the bowtie reads folder.  You'll need several scripts from Harvard's scriptome project: the column_chooser.pl, ID_subtract.pl, 
##  ID_intersect.pl and fasta_picker.pl scripts too.  Also, you'll need the liftOver executable and the droSec1 liftover table installed.

use strict;
use FileHandle;

system 'date';

# This section aligns all the solexa reads to the melanogaster and sechellia genomes
system './bowtie -t -f -v 0 -k 3 -p 4 D_melanogaster_2006 reads/Hyb_oct_r1.fa Hyb_oct_r1_mel.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_melanogaster_2006 reads/Hyb_oct_r2.fa Hyb_oct_r2_mel.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_melanogaster_2006 reads/Mix_oct_r1.fa Mix_oct_r1_mel.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_melanogaster_2006 reads/Mix_oct_r2.fa Mix_oct_r2_mel.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_sechellia_2005 reads/Hyb_oct_r1.fa Hyb_oct_r1_sec.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_sechellia_2005 reads/Hyb_oct_r2.fa Hyb_oct_r2_sec.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_sechellia_2005 reads/Mix_oct_r1.fa Mix_oct_r1_sec.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_sechellia_2005 reads/Mix_oct_r2.fa Mix_oct_r2_sec.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_melanogaster_2006 reads/Hyb_nov_r1.fa Hyb_nov_r1_mel.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_melanogaster_2006 reads/Hyb_nov_r2.fa Hyb_nov_r2_mel.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_melanogaster_2006 reads/Mix_nov_r1.fa Mix_nov_r1_mel.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_melanogaster_2006 reads/Mix_nov_r2.fa Mix_nov_r2_mel.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_sechellia_2005 reads/Hyb_nov_r1.fa Hyb_nov_r1_sec.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_sechellia_2005 reads/Hyb_nov_r2.fa Hyb_nov_r2_sec.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_sechellia_2005 reads/Mix_nov_r1.fa Mix_nov_r1_sec.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_sechellia_2005 reads/Mix_nov_r2.fa Mix_nov_r2_sec.map';


print "Aligned fasta files to Melanogaster and Sechellia genomes\n";

# This section chooses READIDs from the alignment files for future comparison to pick SNPs
system 'perl column_chooser.pl Hyb_oct_r1_mel.map > Hyb_oct_r1_mel.IDs';
system 'perl column_chooser.pl Hyb_oct_r2_mel.map > Hyb_oct_r2_mel.IDs';
system 'perl column_chooser.pl Hyb_oct_r1_sec.map > Hyb_oct_r1_sec.IDs';
system 'perl column_chooser.pl Hyb_oct_r2_sec.map > Hyb_oct_r2_sec.IDs';
system 'perl column_chooser.pl Mix_oct_r1_mel.map > Mix_oct_r1_mel.IDs';
system 'perl column_chooser.pl Mix_oct_r2_mel.map > Mix_oct_r2_mel.IDs';
system 'perl column_chooser.pl Mix_oct_r1_sec.map > Mix_oct_r1_sec.IDs';
system 'perl column_chooser.pl Mix_oct_r2_sec.map > Mix_oct_r2_sec.IDs';
system 'perl column_chooser.pl Hyb_nov_r1_mel.map > Hyb_nov_r1_mel.IDs';
system 'perl column_chooser.pl Hyb_nov_r2_mel.map > Hyb_nov_r2_mel.IDs';
system 'perl column_chooser.pl Hyb_nov_r1_sec.map > Hyb_nov_r1_sec.IDs';
system 'perl column_chooser.pl Hyb_nov_r2_sec.map > Hyb_nov_r2_sec.IDs';
system 'perl column_chooser.pl Mix_nov_r1_mel.map > Mix_nov_r1_mel.IDs';
system 'perl column_chooser.pl Mix_nov_r2_mel.map > Mix_nov_r2_mel.IDs';
system 'perl column_chooser.pl Mix_nov_r1_sec.map > Mix_nov_r1_sec.IDs';
system 'perl column_chooser.pl Mix_nov_r2_sec.map > Mix_nov_r2_sec.IDs';


print "Chose ReadIDs for all bowtie alignment files\n";

# This section compares ReadIDs to pick SNPs, common reads, and SNP-linked reads
system 'perl ID_subtract.pl Hyb_oct_r1_mel.IDs Hyb_oct_r1_sec.IDs > Hyb_oct_r1_mel_SNP.IDs';
system 'perl ID_subtract.pl Hyb_oct_r1_sec.IDs Hyb_oct_r1_mel.IDs > Hyb_oct_r1_sec_SNP.IDs';
system 'perl ID_intersect.pl Hyb_oct_r1_mel.IDs Hyb_oct_r1_sec.IDS > Hyb_oct_r1_common.IDs';
system 'perl ID_subtract.pl Hyb_oct_r2_mel.IDs Hyb_oct_r2_sec.IDs > Hyb_oct_r2_mel_SNP.IDs';
system 'perl ID_subtract.pl Hyb_oct_r2_sec.IDs Hyb_oct_r2_mel.IDs > Hyb_oct_r2_sec_SNP.IDs';
system 'perl ID_intersect.pl Hyb_oct_r2_mel.IDs Hyb_oct_r2_sec.IDS > Hyb_oct_r2_common.IDs';

system 'perl ID_subtract.pl Mix_oct_r1_mel.IDs Mix_oct_r1_sec.IDs > Mix_oct_r1_mel_SNP.IDs';
system 'perl ID_subtract.pl Mix_oct_r1_sec.IDs Mix_oct_r1_mel.IDs > Mix_oct_r1_sec_SNP.IDs';
system 'perl ID_intersect.pl Mix_oct_r1_mel.IDs Mix_oct_r1_sec.IDS > Mix_oct_r1_common.IDs';
system 'perl ID_subtract.pl Mix_oct_r2_mel.IDs Mix_oct_r2_sec.IDs > Mix_oct_r2_mel_SNP.IDs';
system 'perl ID_subtract.pl Mix_oct_r2_sec.IDs Mix_oct_r2_mel.IDs > Mix_oct_r2_sec_SNP.IDs';
system 'perl ID_intersect.pl Mix_oct_r2_mel.IDs Mix_oct_r2_sec.IDS > Mix_oct_r2_common.IDs';

system 'perl ID_intersect.pl Hyb_oct_r1_mel_SNP.IDs Hyb_oct_r2_common.IDS > Hyb_oct_r2_mel_link.IDs';
system 'perl ID_intersect.pl Hyb_oct_r2_mel_SNP.IDs Hyb_oct_r1_common.IDS > Hyb_oct_r1_mel_link.IDs';
system 'perl ID_intersect.pl Hyb_oct_r1_sec_SNP.IDs Hyb_oct_r2_common.IDS > Hyb_oct_r2_sec_link.IDs';
system 'perl ID_intersect.pl Hyb_oct_r2_sec_SNP.IDs Hyb_oct_r1_common.IDS > Hyb_oct_r1_sec_link.IDs';
system 'perl ID_intersect.pl Mix_oct_r1_mel_SNP.IDs Mix_oct_r2_common.IDS > Mix_oct_r2_mel_link.IDs';
system 'perl ID_intersect.pl Mix_oct_r2_mel_SNP.IDs Mix_oct_r1_common.IDS > Mix_oct_r1_mel_link.IDs';
system 'perl ID_intersect.pl Mix_oct_r1_sec_SNP.IDs Mix_oct_r2_common.IDS > Mix_oct_r2_sec_link.IDs';
system 'perl ID_intersect.pl Mix_oct_r2_sec_SNP.IDs Mix_oct_r1_common.IDS > Mix_oct_r1_sec_link.IDs';

system 'perl ID_subtract.pl Hyb_nov_r1_mel.IDs Hyb_nov_r1_sec.IDs > Hyb_nov_r1_mel_SNP.IDs';
system 'perl ID_subtract.pl Hyb_nov_r1_sec.IDs Hyb_nov_r1_mel.IDs > Hyb_nov_r1_sec_SNP.IDs';
system 'perl ID_intersect.pl Hyb_nov_r1_mel.IDs Hyb_nov_r1_sec.IDS > Hyb_nov_r1_common.IDs';
system 'perl ID_subtract.pl Hyb_nov_r2_mel.IDs Hyb_nov_r2_sec.IDs > Hyb_nov_r2_mel_SNP.IDs';
system 'perl ID_subtract.pl Hyb_nov_r2_sec.IDs Hyb_nov_r2_mel.IDs > Hyb_nov_r2_sec_SNP.IDs';
system 'perl ID_intersect.pl Hyb_nov_r2_mel.IDs Hyb_nov_r2_sec.IDS > Hyb_nov_r2_common.IDs';

system 'perl ID_subtract.pl Mix_nov_r1_mel.IDs Mix_nov_r1_sec.IDs > Mix_nov_r1_mel_SNP.IDs';
system 'perl ID_subtract.pl Mix_nov_r1_sec.IDs Mix_nov_r1_mel.IDs > Mix_nov_r1_sec_SNP.IDs';
system 'perl ID_intersect.pl Mix_nov_r1_mel.IDs Mix_nov_r1_sec.IDS > Mix_nov_r1_common.IDs';
system 'perl ID_subtract.pl Mix_nov_r2_mel.IDs Mix_nov_r2_sec.IDs > Mix_nov_r2_mel_SNP.IDs';
system 'perl ID_subtract.pl Mix_nov_r2_sec.IDs Mix_nov_r2_mel.IDs > Mix_nov_r2_sec_SNP.IDs';
system 'perl ID_intersect.pl Mix_nov_r2_mel.IDs Mix_nov_r2_sec.IDS > Mix_nov_r2_common.IDs';

system 'perl ID_intersect.pl Hyb_nov_r1_mel_SNP.IDs Hyb_nov_r2_common.IDS > Hyb_nov_r2_mel_link.IDs';
system 'perl ID_intersect.pl Hyb_nov_r2_mel_SNP.IDs Hyb_nov_r1_common.IDS > Hyb_nov_r1_mel_link.IDs';
system 'perl ID_intersect.pl Hyb_nov_r1_sec_SNP.IDs Hyb_nov_r2_common.IDS > Hyb_nov_r2_sec_link.IDs';
system 'perl ID_intersect.pl Hyb_nov_r2_sec_SNP.IDs Hyb_nov_r1_common.IDS > Hyb_nov_r1_sec_link.IDs';
system 'perl ID_intersect.pl Mix_nov_r1_mel_SNP.IDs Mix_nov_r2_common.IDS > Mix_nov_r2_mel_link.IDs';
system 'perl ID_intersect.pl Mix_nov_r2_mel_SNP.IDs Mix_nov_r1_common.IDS > Mix_nov_r1_mel_link.IDs';
system 'perl ID_intersect.pl Mix_nov_r1_sec_SNP.IDs Mix_nov_r2_common.IDS > Mix_nov_r2_sec_link.IDs';
system 'perl ID_intersect.pl Mix_nov_r2_sec_SNP.IDs Mix_nov_r1_common.IDS > Mix_nov_r1_sec_link.IDs';


print "Compared ReadIDs and chose species-specific(SNP) and SNP-linked Reads\n";

# This section chooses SNP and SNP-linked fasta records for realignment to bowtie (saves time and memory compared to choosing from the original alignment files)

system 'perl fasta_picker.pl Hyb_oct_r1_mel_SNP.IDs reads/Hyb_oct_r1.fa > Hyb_oct_r1_mel_SNP.fa';
system 'perl fasta_picker.pl Hyb_oct_r2_mel_SNP.IDs reads/Hyb_oct_r2.fa > Hyb_oct_r2_mel_SNP.fa';
system 'perl fasta_picker.pl Hyb_oct_r1_mel_link.IDs reads/Hyb_oct_r1.fa > Hyb_oct_r1_mel_link.fa';
system 'perl fasta_picker.pl Hyb_oct_r2_mel_link.IDs reads/Hyb_oct_r2.fa > Hyb_oct_r2_mel_link.fa';
system 'perl fasta_picker.pl Hyb_oct_r1_sec_SNP.IDs reads/Hyb_oct_r1.fa > Hyb_oct_r1_sec_SNP.fa';
system 'perl fasta_picker.pl Hyb_oct_r2_sec_SNP.IDs reads/Hyb_oct_r2.fa > Hyb_oct_r2_sec_SNP.fa';
system 'perl fasta_picker.pl Hyb_oct_r1_sec_link.IDs reads/Hyb_oct_r1.fa > Hyb_oct_r1_sec_link.fa';
system 'perl fasta_picker.pl Hyb_oct_r2_sec_link.IDs reads/Hyb_oct_r2.fa > Hyb_oct_r2_sec_link.fa';
system 'perl fasta_picker.pl Mix_oct_r1_mel_SNP.IDs reads/Mix_oct_r1.fa > Mix_oct_r1_mel_SNP.fa';
system 'perl fasta_picker.pl Mix_oct_r2_mel_SNP.IDs reads/Mix_oct_r2.fa > Mix_oct_r2_mel_SNP.fa';
system 'perl fasta_picker.pl Mix_oct_r1_mel_link.IDs reads/Mix_oct_r1.fa > Mix_oct_r1_mel_link.fa';
system 'perl fasta_picker.pl Mix_oct_r2_mel_link.IDs reads/Mix_oct_r2.fa > Mix_oct_r2_mel_link.fa';
system 'perl fasta_picker.pl Mix_oct_r1_sec_SNP.IDs reads/Mix_oct_r1.fa > Mix_oct_r1_sec_SNP.fa';
system 'perl fasta_picker.pl Mix_oct_r2_sec_SNP.IDs reads/Mix_oct_r2.fa > Mix_oct_r2_sec_SNP.fa';
system 'perl fasta_picker.pl Mix_oct_r1_sec_link.IDs reads/Mix_oct_r1.fa > Mix_oct_r1_sec_link.fa';
system 'perl fasta_picker.pl Mix_oct_r2_sec_link.IDs reads/Mix_oct_r2.fa > Mix_oct_r2_sec_link.fa';

print "Picked SNP and SNP-linked reads from original fasta files\n";

# This section runs bowtie on the SNP and SNP-linked fasta files

system 'mv *.fa reads/';
system './bowtie -t -f -v 0 -k 3 -p 4 D_melanogaster_2006 reads/Mix_oct_r1_mel_link.fa Mix_oct_r1_mel_link.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_melanogaster_2006 reads/Mix_oct_r2_mel_link.fa Mix_oct_r2_mel_link.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_melanogaster_2006 reads/Mix_oct_r1_mel_SNP.fa Mix_oct_r1_mel_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_melanogaster_2006 reads/Mix_oct_r2_mel_SNP.fa Mix_oct_r2_mel_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_sechellia_2005 reads/Mix_oct_r1_sec_link.fa Mix_oct_r1_sec_link.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_sechellia_2005 reads/Mix_oct_r2_sec_link.fa Mix_oct_r2_sec_link.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_sechellia_2005 reads/Mix_oct_r1_sec_SNP.fa Mix_oct_r1_sec_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_sechellia_2005 reads/Mix_oct_r2_sec_SNP.fa Mix_oct_r2_sec_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_melanogaster_2006 reads/Hyb_oct_r1_mel_link.fa Hyb_oct_r1_mel_link.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_melanogaster_2006 reads/Hyb_oct_r2_mel_link.fa Hyb_oct_r2_mel_link.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_melanogaster_2006 reads/Hyb_oct_r1_mel_SNP.fa Hyb_oct_r1_mel_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_melanogaster_2006 reads/Hyb_oct_r2_mel_SNP.fa Hyb_oct_r2_mel_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_sechellia_2005 reads/Hyb_oct_r1_sec_link.fa Hyb_oct_r1_sec_link.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_sechellia_2005 reads/Hyb_oct_r2_sec_link.fa Hyb_oct_r2_sec_link.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_sechellia_2005 reads/Hyb_oct_r1_sec_SNP.fa Hyb_oct_r1_sec_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 4 D_sechellia_2005 reads/Hyb_oct_r2_sec_SNP.fa Hyb_oct_r2_sec_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_melanogaster_2006 reads/Mix_nov_r1_mel_link.fa Mix_nov_r1_mel_link.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_melanogaster_2006 reads/Mix_nov_r2_mel_link.fa Mix_nov_r2_mel_link.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_melanogaster_2006 reads/Mix_nov_r1_mel_SNP.fa Mix_nov_r1_mel_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_melanogaster_2006 reads/Mix_nov_r2_mel_SNP.fa Mix_nov_r2_mel_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_sechellia_2005 reads/Mix_nov_r1_sec_link.fa Mix_nov_r1_sec_link.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_sechellia_2005 reads/Mix_nov_r2_sec_link.fa Mix_nov_r2_sec_link.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_sechellia_2005 reads/Mix_nov_r1_sec_SNP.fa Mix_nov_r1_sec_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_sechellia_2005 reads/Mix_nov_r2_sec_SNP.fa Mix_nov_r2_sec_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_melanogaster_2006 reads/Hyb_nov_r1_mel_link.fa Hyb_nov_r1_mel_link.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_melanogaster_2006 reads/Hyb_nov_r2_mel_link.fa Hyb_nov_r2_mel_link.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_melanogaster_2006 reads/Hyb_nov_r1_mel_SNP.fa Hyb_nov_r1_mel_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_melanogaster_2006 reads/Hyb_nov_r2_mel_SNP.fa Hyb_nov_r2_mel_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_sechellia_2005 reads/Hyb_nov_r1_sec_link.fa Hyb_nov_r1_sec_link.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_sechellia_2005 reads/Hyb_nov_r2_sec_link.fa Hyb_nov_r2_sec_link.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_sechellia_2005 reads/Hyb_nov_r1_sec_SNP.fa Hyb_nov_r1_sec_SNP.map';
system './bowtie -t -f -v 0 -k 3 -p 2 D_sechellia_2005 reads/Hyb_nov_r2_sec_SNP.fa Hyb_nov_r2_sec_SNP.map';

print "Aligned SNP and SNP-linked reads";


print "Aligned SNP and SNP-linked reads\n";

# This section converts the bowtie alignments to BED format files (for liftover and filtering purposes).  Note that the .map files were in 
# 0-base coordinates.  This is the same as BED format, so we add the read length to create the stop coordinate.

my @READS_FILENAME;
$READS_FILENAME[1] = "Hyb_oct_r1_mel_SNP.map";
$READS_FILENAME[2] = "Hyb_oct_r1_mel_link.map";
$READS_FILENAME[3] = "Hyb_oct_r1_sec_SNP.map";
$READS_FILENAME[4] = "Hyb_oct_r1_sec_link.map";
$READS_FILENAME[5] = "Hyb_oct_r2_mel_SNP.map";
$READS_FILENAME[6] = "Hyb_oct_r2_mel_link.map";
$READS_FILENAME[7] = "Hyb_oct_r2_sec_SNP.map";
$READS_FILENAME[8] = "Hyb_oct_r2_sec_link.map";
$READS_FILENAME[9] = "Mix_oct_r1_mel_SNP.map";
$READS_FILENAME[10] = "Mix_oct_r1_mel_link.map";
$READS_FILENAME[11] = "Mix_oct_r1_sec_SNP.map";
$READS_FILENAME[12] = "Mix_oct_r1_sec_link.map";
$READS_FILENAME[13] = "Mix_oct_r2_mel_SNP.map";
$READS_FILENAME[14] = "Mix_oct_r2_mel_link.map";
$READS_FILENAME[15] = "Mix_oct_r2_sec_SNP.map";
$READS_FILENAME[16] = "Mix_oct_r2_sec_link.map";
$READS_FILENAME[17] = "Hyb_nov_r1_mel_SNP.map";
$READS_FILENAME[18] = "Hyb_nov_r1_mel_link.map";
$READS_FILENAME[19] = "Hyb_nov_r1_sec_SNP.map";
$READS_FILENAME[20] = "Hyb_nov_r1_sec_link.map";
$READS_FILENAME[21] = "Hyb_nov_r2_mel_SNP.map";
$READS_FILENAME[22] = "Hyb_nov_r2_mel_link.map";
$READS_FILENAME[23] = "Hyb_nov_r2_sec_SNP.map";
$READS_FILENAME[24] = "Hyb_nov_r2_sec_link.map";
$READS_FILENAME[25] = "Mix_nov_r1_mel_SNP.map";
$READS_FILENAME[26] = "Mix_nov_r1_mel_link.map";
$READS_FILENAME[27] = "Mix_nov_r1_sec_SNP.map";
$READS_FILENAME[28] = "Mix_nov_r1_sec_link.map";
$READS_FILENAME[29] = "Mix_nov_r2_mel_SNP.map";
$READS_FILENAME[30] = "Mix_nov_r2_mel_link.map";
$READS_FILENAME[31] = "Mix_nov_r2_sec_SNP.map";
$READS_FILENAME[32] = "Mix_nov_r2_sec_link.map";



##my $DATA_DIR = "/Users/Joel/Desktop/bowtie-0.9.8";   ####CHANGE_THIS####
####my $LANE = 2;

for(my $ifile=1; $ifile<=32; $ifile++){          #LOOP OVER READ FILES



open(READS, "<$READS_FILENAME[$ifile]") or die("can't open $READS_FILENAME[$ifile]");

my $OUTFILE = $READS_FILENAME[$ifile];
$OUTFILE =~ s/.map/.BED/;

my $BED_OUTPUT = new FileHandle ">$OUTFILE" or die "can't open $OUTFILE";   #reads that hit cons exons


  my $line;
  my $header;
  my @f;
  my $count;
  my $seq;
  my $contig;
  my $startpos;
  my $stoppos;
  my $Strand;
  my $ReadID;

  my $readnum=0;


  while ($line =<READS>) {   
    $readnum++;
    ##printf("$readnum\n");
    chomp($line);
    @f = split(/\s+/,$line);
   	$contig = $f[3];
    $startpos = $f[4];
    $stoppos = $startpos + 37;
    $Strand = $f[2];   
    $ReadID = $f[0];
    $BED_OUTPUT->printf("$contig\t$startpos\t$stoppos\t$ReadID\t".".\t"."$Strand\n");    
  }
printf("Converted $readnum lines from $READS_FILENAME[$ifile] to BED format\n");
}

##  This section lifts the sechellia coordinates to melanogaster space

system './liftOver.MacOSX.ppc Hyb_oct_r1_sec_SNP.BED droSec1ToDm3.over.chain Hyb_oct_r1_sec_SNP_lift.BED Hyb_oct_r1_sec_SNP.unmap';
system './liftOver.MacOSX.ppc Hyb_oct_r1_sec_link.BED droSec1ToDm3.over.chain Hyb_oct_r1_sec_link_lift.BED Hyb_oct_r1_sec_link.unmap';
system './liftOver.MacOSX.ppc Hyb_oct_r2_sec_SNP.BED droSec1ToDm3.over.chain Hyb_oct_r2_sec_SNP_lift.BED Hyb_oct_r2_sec_SNP.unmap';
system './liftOver.MacOSX.ppc Hyb_oct_r2_sec_link.BED droSec1ToDm3.over.chain Hyb_oct_r2_sec_link_lift.BED Hyb_oct_r2_sec_link.unmap';
system './liftOver.MacOSX.ppc Mix_oct_r1_sec_SNP.BED droSec1ToDm3.over.chain Mix_oct_r1_sec_SNP_lift.BED Mix_oct_r1_sec_SNP.unmap';
system './liftOver.MacOSX.ppc Mix_oct_r1_sec_link.BED droSec1ToDm3.over.chain Mix_oct_r1_sec_link_lift.BED Mix_oct_r1_sec_link.unmap';
system './liftOver.MacOSX.ppc Mix_oct_r2_sec_SNP.BED droSec1ToDm3.over.chain Mix_oct_r2_sec_SNP_lift.BED Mix_oct_r2_sec_SNP.unmap';
system './liftOver.MacOSX.ppc Mix_oct_r2_sec_link.BED droSec1ToDm3.over.chain Mix_oct_r2_sec_link_lift.BED Mix_oct_r2_sec_link.unmap';


print "Lifted Sechellia files over to Melanogaster space\n\nFiles are ready for gap-filtering on Galaxy\n";
system 'date';

