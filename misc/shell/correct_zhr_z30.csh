#!/bin/csh

#mapping bias

perl correctMapBias_cut20.pl zhr z30 ../mel_mel_data/mRNA-Seq/zhr/zhr.mosaik.exons.txt ../mel_mel_data/mRNA-Seq/z30/z30.mosaik.exons.txt ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mosaik.exons.txt ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mosaik.exons.txt ../mel_mel_data mosaik.exons.txt &

perl correctMapBias_cut20.pl zhr z30 ../mel_mel_data/mRNA-Seq/zhr/zhr.mosaik.genes.txt ../mel_mel_data/mRNA-Seq/z30/z30.mosaik.genes.txt ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mosaik.genes.txt ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mosaik.genes.txt ../mel_mel_data mosaik.genes.txt &

perl correctMapBias_cut20.pl zhr z30 ../mel_mel_data/mRNA-Seq/zhr/zhr_v2.mosaik.exons.txt ../mel_mel_data/mRNA-Seq/z30/z30_v2.mosaik.exons.txt ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_v2.mosaik.exons.txt ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_v2.mosaik.exons.txt ../mel_mel_data/mRNA-Seq mosaik.exons.txt &

#factorParS1 = 1
#factorParS2 = 0.596351976946768
#factorHyb1S1 = 0.970766791705004
#factorHyb1S2 = 1
#factorHyb2S1 = 0.973897657991578
#factorHyb2S2 = 1

perl correctMapBias_cut20.pl zhr z30 ../mel_mel_data/mRNA-Seq/zhr/zhr_v2.mosaik.genes.txt ../mel_mel_data/mRNA-Seq/z30/z30_v2.mosaik.genes.txt ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_v2.mosaik.genes.txt ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_v2.mosaik.genes.txt ../mel_mel_data/mRNA-Seq mosaik.genes.txt &

#factorParS1 = 1
#factorParS2 = 0.611867957501099
#factorHyb1S1 = 0.970560306204448
#factorHyb1S2 = 1
#factorHyb2S1 = 0.973918764627187
#factorHyb2S2 = 1


#sequencing depth

#perl correctSeqDepth_cut20.pl zhr z30 ../mel_mel_data/mRNA-Seq/zhr/zhr.mosaik.exons.txt 16464075 ../mel_mel_data/mRNA-Seq/z30/z30.mosaik.exons.txt 21806797 ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mosaik.exons.txt 31432754 ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mosaik.exons.txt 31439998 ../mel_mel_data mosaik.exons.txt &

#perl correctSeqDepth_cut20.pl zhr z30 ../mel_mel_data/mRNA-Seq/zhr/zhr.mosaik.genes.txt 16464075 ../mel_mel_data/mRNA-Seq/z30/z30.mosaik.genes.txt 21806797 ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mosaik.genes.txt 31432754 ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mosaik.genes.txt 31439998 ../mel_mel_data mosaik.genes.txt &


#perl correctSD_cut20_combine_sepFile2.pl zhr z30 zhr_mRNA ../mel_mel_data/mRNA-Seq/zhr/zhr_v2.mosaik.genes.txt 16464075 z30_mRNA ../mel_mel_data/mRNA-Seq/z30/z30_v2.mosaik.genes.txt 21806797 ../mel_mel_data/mRNA-Seq/zhr+z30_mRNA_SD_cut20.mosaik.genes.txt

#perl correctSD_cut20_combine_sepFile2.pl zhr z30 zhr_mRNA ../mel_mel_data/mRNA-Seq/zhr/zhr_v2.mosaik.exons.txt 16464075 z30_mRNA ../mel_mel_data/mRNA-Seq/z30/z30_v2.mosaik.exons.txt 21806797 ../mel_mel_data/mRNA-Seq/zhr+z30_mRNA_SD_cut20.mosaik.exons.txt


perl combine_files.pl zhr z30 zhrXz30_mRNA ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_v3.mosaik.genes.txt z30Xzhr_mRNA ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_v3.mosaik.genes.txt ../mel_mel_data/mRNA-Seq/zhrXz30_z30Xzhr_mRNA.mosaik.genes.txt

perl combine_files.pl zhr z30 zhrXz30_mRNA ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_v3.mosaik.exons.txt z30Xzhr_mRNA ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_v3.mosaik.exons.txt ../mel_mel_data/mRNA-Seq/zhrXz30_z30Xzhr_mRNA.mosaik.exons.txt


perl combine_files.pl zhr z30 zhrXz30_mRNA ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_v3.mosaik.genes.txt z30Xzhr_mRNA ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_v3.mosaik.genes.txt ../mel_mel_data/mRNA-Seq/zhrXz30_z30Xzhr_mRNA.mosaik.genes.txt

perl combine_files.pl zhr z30 zhrXz30_mRNA ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_v3.mosaik.exons.txt z30Xzhr_mRNA ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_v3.mosaik.exons.txt ../mel_mel_data/mRNA-Seq/zhrXz30_z30Xzhr_mRNA.mosaik.exons.txt


#sequencing depth: gDNA

perl correctSD_combine_sepFile2_gDNA.pl zhr z30 zhr_gDNA ../mel_mel_data/Resequencing/zhr/zhr_v2.mosaik.genes.txt 42743562 z30_gDNA ../mel_mel_data/Resequencing/z30/z30_v2.mosaik.genes.txt 25863911 ../mel_mel_data/Resequencing/zhr_z30_gDNA_SD.mosaik.genes.txt

perl correctSD_combine_sepFile2_gDNA.pl zhr z30 zhr_gDNA ../mel_mel_data/Resequencing/zhr/zhr_v2.mosaik.exons.txt 42743562 z30_gDNA ../mel_mel_data/Resequencing/z30/z30_v2.mosaik.exons.txt 25863911 ../mel_mel_data/Resequencing/zhr_z30_gDNA_SD.mosaik.exons.txt




# zhrXz30 and z30Xzhr

perl combine_files.pl zhr z30 zhrXz30_mRNA ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_v2.mosaik.exons.txt z30Xzhr_mRNA ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_v2.mosaik.exons.txt ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/no_gap_filtering/zhrXz30_z30Xzhr_mRNA.mosaik.exons.txt

perl combine_files.pl zhr z30 zhrXz30_mRNA ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_v2.mosaik.genes.txt z30Xzhr_mRNA ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_v2.mosaik.genes.txt ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/no_gap_filtering/zhrXz30_z30Xzhr_mRNA.mosaik.genes.txt

# zhr+z30 and zhrXz30

perl combine_files.pl zhr z30 zhr+z30_mRNA ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/zhr+z30/zhr+z30_mRNA_SD.mosaik.exons.txt zhrXz30_mRNA ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_v2.mosaik.exons.txt ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/no_gap_filtering/zhr+z30_zhrXz30_mRNA.mosaik.exons.txt

perl combine_files.pl zhr z30 zhr+z30_mRNA ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/zhr+z30/zhr+z30_mRNA_SD.mosaik.genes.txt zhrXz30_mRNA ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_v2.mosaik.genes.txt ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/no_gap_filtering/zhr+z30_zhrXz30_mRNA.mosaik.genes.txt

# zhr+z30 and z30Xzhr

perl combine_files.pl zhr z30 zhr+z30_mRNA ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/zhr+z30/zhr+z30_mRNA_SD.mosaik.exons.txt z30Xzhr_mRNA ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_v2.mosaik.exons.txt ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/no_gap_filtering/zhr+z30_z30Xzhr_mRNA.mosaik.exons.txt

perl combine_files.pl zhr z30 zhr+z30_mRNA ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/zhr+z30/zhr+z30_mRNA_SD.mosaik.genes.txt z30Xzhr_mRNA ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_v2.mosaik.genes.txt ../../../../Volumes/Wittkopp_data/02012012_regEvol_dataFreeze/mel_mel_data/mRNA-Seq/no_gap_filtering/zhr+z30_z30Xzhr_mRNA.mosaik.genes.txt
