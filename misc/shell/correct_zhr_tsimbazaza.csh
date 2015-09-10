#!/bin/csh

#mapping bias

perl correctMapBias_cut20.0.2.pl zhr tsimbazaza ../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mosaik.exons.txt ../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mosaik.exons.txt ../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mosaik.exons.txt ../mel_sim_data mosaik.exons.txt 

perl correctMapBias_cut20.0.2.pl zhr tsimbazaza ../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mosaik.genes.txt ../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mosaik.genes.txt ../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mosaik.genes.txt ../mel_sim_data mosaik.genes.txt


perl correctMB_cut20_sepPar.pl zhr tsimbazaza ../mel_sim_data/mRNA-Seq/zhr/zhr.mosaik.exons.txt ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mosaik.exons.txt ../mel_sim_data mosaik.exons.txt

perl correctMB_cut20_sepPar.pl zhr tsimbazaza ../mel_sim_data/mRNA-Seq/zhr/zhr.mosaik.genes.txt ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mosaik.genes.txt ../mel_sim_data mosaik.genes.txt

#sequencing depth

perl correctSeqDepth_cut20.0.2.pl zhr tsimbazaza ../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mosaik.exons.txt ../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mosaik.exons.txt 23929242 ../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mosaik.exons.txt 25875801 ../mel_sim_data mosaik.exons.txt 

perl correctSeqDepth_cut20.0.2.pl zhr tsimbazaza ../mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mosaik.genes.txt ../mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mosaik.genes.txt 23929242 ../mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mosaik.genes.txt 25875801 ../mel_sim_data mosaik.genes.txt 


perl correctSD_cut20_sepPar.pl zhr tsimbazaza ../mel_sim_data/mRNA-Seq/zhr/zhr.mosaik.exons.txt 16464075 ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mosaik.exons.txt 18006673 ../mel_sim_data mosaik.exons.txt

perl correctSD_cut20_sepPar.pl zhr tsimbazaza ../mel_sim_data/mRNA-Seq/zhr/zhr.mosaik.genes.txt 16464075 ../mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza.mosaik.genes.txt 18006673 ../mel_sim_data mosaik.genes.txt
