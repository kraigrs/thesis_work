#!/bin/sh

# lift differentiating sites in regulatory regions

# sim
 
liftOver ../../validation/dm3_candidate_regs.SNPs.bed ../../genomes/dm3/dm3ToDroSim1.over.chain ../../validation/droSim1_candidate_regs.SNPs.bed ../../validation/droSim1_candidate_regs.SNPs.un

liftOver ../../validation/droSim1_candidate_regs.SNPs.bed ../../genomes/tsimbazaza/droSim1ToDroSim1_snp_indel.over.chain ../../validation/sim_candidate_regs.SNPs.bed ../../validation/sim_candidate_regs.SNPs.un

liftOver ../../validation/droSim1_candidate_regs.SNPs.un ../../genomes/dm3/dm3ToSim_extra.over.chain ../../validation/sim_extra_candidate_regs.SNPs.bed ../../validation/sim_extra_candidate_regs.SNPs.un

# sec
 
liftOver ../../validation/dm3_candidate_regs.SNPs.bed ../../genomes/dm3/dm3ToDroSec1.over.chain ../../validation/droSec1_candidate_regs.SNPs.bed ../../validation/droSec1_candidate_regs.SNPs.un

liftOver ../../validation/droSec1_candidate_regs.SNPs.bed ../../genomes/droSec1/droSec1ToDroSec1_snp_indel.over.chain ../../validation/sec_candidate_regs.SNPs.bed ../../validation/sec_candidate_regs.SNPs.un

liftOver ../../validation/droSec1_candidate_regs.SNPs.un ../../genomes/dm3/dm3ToSec_extra.over.chain ../../validation/sec_extra_candidate_regs.SNPs.bed ../../validation/sec_extra_candidate_regs.SNPs.un
