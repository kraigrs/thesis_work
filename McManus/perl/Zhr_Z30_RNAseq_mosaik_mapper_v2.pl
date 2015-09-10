#!/usr/bin/perl
##  Joel McManus 07/23/2010
##  This script aligns RNA-seq data to find snp genomic and junction hits from Zhr and Z30.
###  !!!! NOTES !!!! Check the readIDs to make sure same lane, different date doesn't have same ID naming scheme!!!!!

## for Yang Bing: 
## I've edited this script to remove the splice junction mapping and trans-spliced junction mapping

use strict;
use FileHandle;
system "date";
printf("\n**********\n\tAligning reads to genome and splice junctions\n**********\n\n");
# This section aligns all the RNA-seq reads to the melanogaster zhr and z30 genomes and junctions. Reads that don't map to genome are written into "nomap" files and then aligned to junctions. Reads that don't map to genome or junction are then aligned to trans-spliced junctions.
my @readfiles = ("3144CAAXX_l6_r1_zhr_seq.txt", "3144CAAXX_l6_r2_zhr_seq.txt", "61FOVAAXX_l4_r2_z30_seq.txt", "61FOVAAXX_l4_r1_z30_seq.txt", "61F0CAAXX_l3_r2_z30xzhr_seq.txt", "61F0CAAXX_l3_r1_z30xzhr_seq.txt", "61F0CAAXX_l2_r2_zhrxz30_seq.txt", "61F0CAAXX_l2_r1_zhrxz30_seq.txt");
foreach (@readfiles) {
	my $nomap1 = $_; $nomap1 =~ s/seq.txt/nomap1.fq/;
	my $nomap2 = $_; $nomap2 =~ s/seq.txt/nomap2.fq/;
	my $nomap3 = $_; $nomap3 =~ s/seq.txt/nomap3.fq/;
	my $nomap4 = $_; $nomap4 =~ s/seq.txt/nomap4.fq/;

	my $map_zhr = $_; $map_zhr =~ s/seq.txt/zhr.dat/;
	my $map_z30 = $_; $map_z30 =~ s/seq.txt/z30.dat/;
	my $snp_zhr = $_; $snp_zhr =~ s/seq.txt/zhr_snp.dat/;
	my $snp_z30 = $_; $snp_z30 =~ s/seq.txt/z30_snp.dat/;

	my $trim = 76;
	my $mosaikreads;
	$mosaikreads = $_;
	printf("$_ and $mosaikreads\n\n");
	$mosaikreads =~ s/seq.txt/input.dat/;
	my $startfile;
	# Align and find snp reads to genomic and junction databases.
	while ($trim >= 37) {
		if ($trim == 76) {
			$startfile = $_;
		} else {
			$startfile = "$_"."_$trim";
		}
		printf("\n*********\nAligning reads in file $startfile to D. mel zhr and D. mel z30 genome databases.\n*********\n");

		system "MosaikBuild -q $startfile -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $map_zhr"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_AMB.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/zhr_AMB_13 -mhp 100 -act 20 -rur $nomap1";

		system "MosaikBuild -q $nomap3 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snp_z30"."_$trim"." -ia ~/mosaik-aligner/targets/z30_AMB.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/z30_AMB_13 -mhp 100 -act 20 -rur $nomap1";

		system "MosaikBuild -q $startfile -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $map_z30"."_$trim"." -ia ~/mosaik-aligner/targets/z30_AMB.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/z30_AMB_13 -mhp 100 -act 20 -rur $nomap1";

		system "MosaikBuild -q $nomap3 -out $mosaikreads -st illumina";
		system "MosaikAligner -in $mosaikreads -out $snp_zhr"."_$trim"." -ia ~/mosaik-aligner/targets/zhr_AMB.dat -hs 13 -mm 0 -p 24 -j  ~/mosaik-aligner/targets/zhr_AMB_13 -mhp 100 -act 20 -rur $nomap1";


		$trim -= 13;
		if ($trim > 24) {
			printf("\n*********\nTrimming non-hitter reads to length=$trim for remapping.\n\n");
			system "perl fastq_trimmer.pl $nomap2 $trim > "."$_"."_$trim";	## trim reads for next round.
			}
		} ### Repeat above with trimmed reads and junctions.
} # align next reads file
system "date";
