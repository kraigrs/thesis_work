#!/bin/csh

perl correctMapBias_cut20.0.2.pl dm3 droSec1 ../mel_sec_data/mRNA-Seq/Mix.mRNA.mosaik.exonout.txt ../mel_mel_data/mRNA-Seq/z30/z30.mRNA.mosaik.exonout.txt ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mRNA.mosaik.exonout.txt ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mRNA.mosaik.exonout.txt ../mel_sec_data mosaik.exonout.txt

perl correctMapBias_cut20.0.2.pl dm3 droSec1 ../mel_sec_data/mRNA-Seq/Mix.mRNA.mosaik.geneout.txt ../mel_mel_data/mRNA-Seq/z30/z30.mRNA.mosaik.geneout.txt ../mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mRNA.mosaik.exonout.txt ../mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mRNA.mosaik.geneout.txt ../mel_sec_data mosaik.geneout.txt
