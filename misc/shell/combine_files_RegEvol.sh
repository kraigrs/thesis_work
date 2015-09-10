#!/bin/sh

# mel_mel

# zhr and z30 separate parents

perl ../perlScripts/combine_files.pl zhr z30 zhr_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhr/zhr_mRNA.mosaik.genes.txt z30_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/z30/z30_mRNA.mosaik.genes.txt ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhr_z30_mRNA.mosaik.genes.txt

perl ../perlScripts/combine_files.pl zhr z30 zhr_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhr/zhr_mRNA.mosaik.exons.txt z30_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/z30/z30_mRNA.mosaik.exons.txt ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhr_z30_mRNA.mosaik.exons.txt

# zhr+z30 and zhrXz30

perl ../perlScripts/combine_files.pl zhr z30 zhr+z30_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhr+z30/zhr+z30_mRNA_SD.mosaik.genes.txt zhrXz30_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_mRNA.mosaik.genes.txt ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhr+z30_zhrXz30_mRNA.mosaik.genes.txt

perl ../perlScripts/combine_files.pl zhr z30 zhr+z30_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhr+z30/zhr+z30_mRNA_SD.mosaik.exons.txt zhrXz30_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_mRNA.mosaik.exons.txt ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhr+z30_zhrXz30_mRNA.mosaik.exons.txt

# zhr+z30 and z30Xzhr

perl ../perlScripts/combine_files.pl zhr z30 zhr+z30_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhr+z30/zhr+z30_mRNA_SD.mosaik.genes.txt z30Xzhr_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_mRNA.mosaik.genes.txt ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhr+z30_z30Xzhr_mRNA.mosaik.genes.txt

perl ../perlScripts/combine_files.pl zhr z30 zhr+z30_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhr+z30/zhr+z30_mRNA_SD.mosaik.exons.txt z30Xzhr_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_mRNA.mosaik.exons.txt ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhr+z30_z30Xzhr_mRNA.mosaik.exons.txt

# zhrXz30 and z30Xzhr

perl ../perlScripts/combine_files.pl zhr z30 zhrXz30_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_mRNA.mosaik.genes.txt z30Xzhr_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_mRNA.mosaik.genes.txt ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhrXz30_z30Xzhr_mRNA.mosaik.genes.txt

perl ../perlScripts/combine_files.pl zhr z30 zhrXz30_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30_mRNA.mosaik.exons.txt z30Xzhr_mRNA ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr_mRNA.mosaik.exons.txt ../../Desktop/data_freeze_12112012/mel_mel_data/mRNA-Seq/zhrXz30_z30Xzhr_mRNA.mosaik.exons.txt

# gDNA

perl ../perlScripts/combine_files.pl zhr z30 zhr_gDNA ../../Desktop/data_freeze_12112012/mel_mel_data/Resequencing/zhr/zhr_gDNA.mosaik.genes.txt z30_gDNA ../../Desktop/data_freeze_12112012/mel_mel_data/Resequencing/z30/z30_gDNA.mosaik.genes.txt ../../Desktop/data_freeze_12112012/mel_mel_data/Resequencing/zhr_z30_gDNA.mosaik.genes.txt

perl ../perlScripts/combine_files.pl zhr z30 zhr_gDNA ../../Desktop/data_freeze_12112012/mel_mel_data/Resequencing/zhr/zhr_gDNA.mosaik.exons.txt z30_gDNA ../../Desktop/data_freeze_12112012/mel_mel_data/Resequencing/z30/z30_gDNA.mosaik.exons.txt ../../Desktop/data_freeze_12112012/mel_mel_data/Resequencing/zhr_z30_gDNA.mosaik.exons.txt

# mel_sim

# zhr and tsimbazaza separate parents

perl ../perlScripts/combine_files.pl zhr tsimbazaza zhr_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhr/zhr_mRNA.mosaik.genes.txt tsimbazaza_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza_mRNA.mosaik.genes.txt ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhr_tsimbazaza_mRNA.mosaik.genes.txt

perl ../perlScripts/combine_files.pl zhr tsimbazaza zhr_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhr/zhr_mRNA.mosaik.exons.txt tsimbazaza_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/tsimbazaza/tsimbazaza_mRNA.mosaik.exons.txt ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhr_tsimbazaza_mRNA.mosaik.exons.txt

# zhr+tsimbazaza and zhrXtsimbazaza

perl ../perlScripts/combine_files.pl zhr tsimbazaza zhr+tsimbazaza_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza_mRNA.mosaik.genes.txt zhrXtsimbazaza_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza_mRNA.mosaik.genes.txt ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhr+tsimbazaza_zhrXtsimbazaza_mRNA.mosaik.genes.txt

perl ../perlScripts/combine_files.pl zhr tsimbazaza zhr+tsimbazaza_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza_mRNA.mosaik.exons.txt zhrXtsimbazaza_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza_mRNA.mosaik.exons.txt ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhr+tsimbazaza_zhrXtsimbazaza_mRNA.mosaik.exons.txt

# zhr+tsimbazaza and tsimbazazaXzhr

perl ../perlScripts/combine_files.pl zhr tsimbazaza zhr+tsimbazaza_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza_mRNA.mosaik.genes.txt tsimbazazaXzhr_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr_mRNA.mosaik.genes.txt ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhr+tsimbazaza_tsimbazazaXzhr_mRNA.mosaik.genes.txt

perl ../perlScripts/combine_files.pl zhr tsimbazaza zhr+tsimbazaza_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza_mRNA.mosaik.exons.txt tsimbazazaXzhr_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr_mRNA.mosaik.exons.txt ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhr+tsimbazaza_tsimbazazaXzhr_mRNA.mosaik.exons.txt

# zhrXtsimbazaza and tsimbazazaXzhr

perl ../perlScripts/combine_files.pl zhr tsimbazaza zhrXtsimbazaza_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza_mRNA.mosaik.genes.txt tsimbazazaXzhr_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr_mRNA.mosaik.genes.txt ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhrXtsimbazaza_tsimbazazaXzhr_mRNA.mosaik.genes.txt

perl ../perlScripts/combine_files.pl zhr tsimbazaza zhrXtsimbazaza_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza_mRNA.mosaik.exons.txt tsimbazazaXzhr_mRNA ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr_mRNA.mosaik.exons.txt ../../Desktop/data_freeze_12112012/mel_sim_data/mRNA-Seq/zhrXtsimbazaza_tsimbazazaXzhr_mRNA.mosaik.exons.txt

# gDNA

perl ../perlScripts/combine_files.pl zhr tsimbazaza zhr_gDNA ../../Desktop/data_freeze_12112012/mel_sim_data/Resequencing/zhr/zhr_gDNA.mosaik.genes.txt tsimbazaza_gDNA ../../Desktop/data_freeze_12112012/mel_sim_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mosaik.genes.txt ../../Desktop/data_freeze_12112012/mel_sim_data/Resequencing/zhr_tsimbazaza_gDNA.mosaik.genes.txt

perl ../perlScripts/combine_files.pl zhr tsimbazaza zhr_gDNA ../../Desktop/data_freeze_12112012/mel_sim_data/Resequencing/zhr/zhr_gDNA.mosaik.exons.txt tsimbazaza_gDNA ../../Desktop/data_freeze_12112012/mel_sim_data/Resequencing/tsimbazaza/tsimbazaza_gDNA.mosaik.exons.txt ../../Desktop/data_freeze_12112012/mel_sim_data/Resequencing/zhr_tsimbazaza_gDNA.mosaik.exons.txt

# sim_sec

# tsimbazaza and droSec1 separate parents

#perl ../perlScripts/combine_files.pl tsimbazaza droSec1 tsimbazaza_mRNA ../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza_v2.mosaik.genes.txt droSec1_mRNA ../sim_sec_data/mRNA-Seq/droSec1/droSec1_v2.mosaik.genes.txt ../sim_sec_data/mRNA-Seq/tsimbazaza_droSec1_mRNA.mosaik.genes.txt

#perl ../perlScripts/combine_files.pl tsimbazaza droSec1 tsimbazaza_mRNA ../sim_sec_data/mRNA-Seq/tsimbazaza/tsimbazaza_v2.mosaik.exons.txt droSec1_mRNA ../sim_sec_data/mRNA-Seq/droSec1/droSec1_v2.mosaik.exons.txt ../sim_sec_data/mRNA-Seq/tsimbazaza_droSec1_mRNA.mosaik.exons.txt

# tsimbazaza+droSec1 and tsimbazazaXdroSec1

#perl ../perlScripts/combine_files.pl tsimbazaza droSec1 tsimbazaza+droSec1_mRNA ../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1_v2.mosaik.genes.txt tsimbazazaXdroSec1_mRNA ../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1_v2.mosaik.genes.txt ../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1_tsimbazazaXdroSec1_mRNA.mosaik.genes.txt

#perl ../perlScripts/combine_files.pl tsimbazaza droSec1 tsimbazaza+droSec1_mRNA ../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1_v2.mosaik.exons.txt tsimbazazaXdroSec1_mRNA ../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1_v2.mosaik.exons.txt ../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1_tsimbazazaXdroSec1_mRNA.mosaik.exons.txt

# tsimbazaza+droSec1 and droSec1Xtsimbazaza

#perl ../perlScripts/combine_files.pl tsimbazaza droSec1 tsimbazaza+droSec1_mRNA ../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1_v2.mosaik.genes.txt droSec1Xtsimbazaza_mRNA ../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza_v2.mosaik.genes.txt ../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1_droSec1Xtsimbazaza_mRNA.mosaik.genes.txt

#perl ../perlScripts/combine_files.pl tsimbazaza droSec1 tsimbazaza+droSec1_mRNA ../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1/tsimbazaza+droSec1_v2.mosaik.exons.txt droSec1Xtsimbazaza_mRNA ../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza_v2.mosaik.exons.txt ../sim_sec_data/mRNA-Seq/tsimbazaza+droSec1_droSec1Xtsimbazaza_mRNA.mosaik.exons.txt

# tsimbazazaXdroSec1 and droSec1Xtsimbazaza

#perl ../perlScripts/combine_files.pl tsimbazaza droSec1 tsimbazazaXdroSec1_mRNA ../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1_v2.mosaik.genes.txt droSec1Xtsimbazaza_mRNA ../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza_v2.mosaik.genes.txt ../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1_droSec1Xtsimbazaza_mRNA.mosaik.genes.txt

#perl ../perlScripts/combine_files.pl tsimbazaza droSec1 tsimbazazaXdroSec1_mRNA ../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1/tsimbazazaXdroSec1_v2.mosaik.exons.txt droSec1Xtsimbazaza_mRNA ../sim_sec_data/mRNA-Seq/droSec1Xtsimbazaza/droSec1Xtsimbazaza_v2.mosaik.exons.txt ../sim_sec_data/mRNA-Seq/tsimbazazaXdroSec1_droSec1Xtsimbazaza_mRNA.mosaik.exons.txt

# gDNA

#perl ../perlScripts/combine_files.pl tsimbazaza droSec1 tsimbazaza_gDNA ../sim_sec_data/Resequencing/tsimbazaza/tsimbazaza_gDNA_v2.mosaik.genes.txt droSec1_gDNA ../sim_sec_data/Resequencing/droSec1/droSec1_gDNA_v2.mosaik.genes.txt  ../sim_sec_data/Resequencing/tsimbazaza_droSec1_gDNA.mosaik.genes.txt

#perl ../perlScripts/combine_files.pl tsimbazaza droSec1 tsimbazaza_gDNA ../sim_sec_data/Resequencing/tsimbazaza/tsimbazaza_gDNA_v2.mosaik.exons.txt droSec1_gDNA ../sim_sec_data/Resequencing/droSec1/droSec1_gDNA_v2.mosaik.exons.txt  ../sim_sec_data/Resequencing/tsimbazaza_droSec1_gDNA.mosaik.exons.txt