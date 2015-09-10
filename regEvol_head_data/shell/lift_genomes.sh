#!/bin/sh

var1=../../genomes/lifted/constitutive_regions_overlap_filtered2 var2=../../genomes/dm3/dm3ToZhr_AMBIG_all.chain var3=zhr qsub liftOver.sh

var1=../../genomes/lifted/constitutive_regions_overlap_filtered2 var2=../../genomes/dm3/dm3ToZ30_AMBIG_all.chain var3=z30 qsub liftOver.sh

var1=../../genomes/lifted/constitutive_regions_overlap_filtered2 var2=../../genomes/dm3/dm3ToDroSim1.over.chain var3=../../genomes/tsimbazaza/droSim1ToDroSim1_snp_indel.over.chain var4=../../genomes/dm3/dm3ToSim_extra.over.chain var5=sim qsub liftOver2.sh

var1=../../genomes/lifted/constitutive_regions_overlap_filtered2 var2=../../genomes/dm3/dm3ToDroSec1.over.chain var3=../../genomes/droSec1/droSec1ToDroSec1_snp_indel.over.chain var4=../../genomes/dm3/dm3ToSec_extra.over.chain var5=sec qsub liftOver2.sh