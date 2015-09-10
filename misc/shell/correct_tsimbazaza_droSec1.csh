#!/bin/csh

#mapping bias

perl correctMapBias_cut20.0.2.pl tsimbazaza droSec1 ../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1.mRNA.mosaik.exons.txt ../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1.mRNA.mosaik.exons.txt ../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza.mRNA.mosaik.exons.txt ../sim_sec_data mosaik.exons.txt 

perl correctMapBias_cut20.0.2.pl tsimbazaza droSec1 ../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1.mRNA.mosaik.genes.txt ../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1.mRNA.mosaik.genes.txt ../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza.mRNA.mosaik.genes.txt ../sim_sec_data mosaik.genes.txt 


perl correctMB_cut20_sepPar2.pl tsimbazaza droSec1 ../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mosaik.exons.txt ../sim_sec_data/mRNA-Seq/droSec1/droSec1.mosaik.exons.txt ../sim_sec_data mosaik.exons.txt

perl correctMB_cut20_sepPar2.pl tsimbazaza droSec1 ../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mosaik.genes.txt ../sim_sec_data/mRNA-Seq/droSec1/droSec1.mosaik.genes.txt ../sim_sec_data mosaik.genes.txt

# using only parent-specific alignment
perl correctMB_cut20_sepPar.pl tsimbazaza droSec1 ../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mosaik.exons.sepPar.txt ../sim_sec_data/mRNA-Seq/droSec1/droSec1.mosaik.exons.sepPar.txt ../sim_sec_data mosaik.exons.txt

perl correctMB_cut20_sepPar.pl tsimbazaza droSec1 ../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mosaik.genes.sepPar.txt ../sim_sec_data/mRNA-Seq/droSec1/droSec1.mosaik.genes.sepPar.txt ../sim_sec_data mosaik.genes.txt

#sequencing depth

perl correctSeqDepth_cut20.0.2.pl tsimbazaza droSec1 ../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1.mRNA.mosaik.exons.txt ../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1.mRNA.mosaik.exons.txt 19787136 ../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza.mRNA.mosaik.exons.txt 20059660 ../sim_sec_data mosaik.exons.txt 

perl correctSeqDepth_cut20.0.2.pl tsimbazaza droSec1 ../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1.mRNA.mosaik.genes.txt ../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1.mRNA.mosaik.genes.txt 19787136 ../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza.mRNA.mosaik.genes.txt 20059660 ../sim_sec_data mosaik.genes.txt 


perl correctSD_cut20_sepPar2.pl tsimbazaza droSec1 ../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mosaik.exons.sepPar.txt 18006673 ../sim_sec_data/mRNA-Seq/droSec1/droSec1.mosaik.exons.sepPar.txt 15817452 ../sim_sec_data mosaik.exons.txt

perl correctSD_cut20_sepPar2.pl tsimbazaza droSec1 ../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mosaik.genes.sepPar.txt 18006673 ../sim_sec_data/mRNA-Seq/droSec1/droSec1.mosaik.genes.sepPar.txt 15817452 ../sim_sec_data mosaik.genes.txt

# using only parent-specific alignment
perl correctSD_cut20_sepPar.pl tsimbazaza droSec1 ../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mosaik.exons.sepPar.txt 18006673 ../sim_sec_data/mRNA-Seq/droSec1/droSec1.mosaik.exons.sepPar.txt 15817452 ../sim_sec_data mosaik.exons.txt

perl correctSD_cut20_sepPar.pl tsimbazaza droSec1 ../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza.mosaik.genes.sepPar.txt 18006673 ../sim_sec_data/mRNA-Seq/droSec1/droSec1.mosaik.genes.sepPar.txt 15817452 ../sim_sec_data mosaik.genes.txt
