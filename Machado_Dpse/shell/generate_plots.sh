#!/bin/sh

# regulatory divergence

R --vanilla --args ../../data/TL_Toro1_H6_F_carcass.genes_down_ASE.txt ../../plots TL_Toro1_H6_F_carcass < ../R/regulatory_divergence.R > ../R/TL_Toro1_H6_F_carcass_reg_div.out

R --vanilla --args ../../data/TL_Toro1_H5_F_carcass.genes_down_ASE.txt ../../plots TL_Toro1_H5_F_carcass < ../R/regulatory_divergence.R > ../R/TL_Toro1_H5_F_carcass_reg_div.out

R --vanilla --args ../../data/TL_Toro1_H6_M_carcass.genes_down_ASE.txt ../../plots TL_Toro1_H6_M_carcass < ../R/regulatory_divergence_no_X.R > ../R/TL_Toro1_H6_M_carcass_reg_div.out

R --vanilla --args ../../data/TL_Toro1_H5_M_carcass.genes_down_ASE.txt ../../plots TL_Toro1_H5_M_carcass < ../R/regulatory_divergence_no_X.R > ../R/TL_Toro1_H5_M_carcass_reg_div.out

R --vanilla --args ../../data/TL_Toro1_H6_ovaries.genes_down_ASE.txt ../../plots TL_Toro1_H6_ovaries < ../R/regulatory_divergence.R > ../R/TL_Toro1_H6_ovaries_reg_div.out

R --vanilla --args ../../data/TL_Toro1_H5_ovaries.genes_down_ASE.txt ../../plots TL_Toro1_H5_ovaries < ../R/regulatory_divergence.R > ../R/TL_Toro1_H5_ovaries_reg_div.out

R --vanilla --args ../../data/TL_Toro1_H6_testes.genes_down_ASE.txt ../../plots TL_Toro1_H6_testes < ../R/regulatory_divergence_no_X.R > ../R/TL_Toro1_H6_testes_reg_div.out

R --vanilla --args ../../data/TL_Toro1_H5_testes.genes_down_ASE.txt ../../plots TL_Toro1_H5_testes < ../R/regulatory_divergence_no_X.R > ../R/TL_Toro1_H5_testes_reg_div.out

# mode of inheritance

R --vanilla --args ../../data/TL_F_carcass.genes_down_total.txt ../../data/Toro1_F_carcass.genes_down_total.txt ../../data/H6_F_carcass.genes_down_total.txt ../../plots TL_Toro1_H6_F_carcass < ../R/mode_of_inheritance.R > ../R/TL_Toro1_H6_F_carcass_MOI.out

R --vanilla --args ../../data/TL_F_carcass.genes_down_total.txt ../../data/Toro1_F_carcass.genes_down_total.txt ../../data/H5_F_carcass.genes_down_total.txt ../../plots TL_Toro1_H5_F_carcass < ../R/mode_of_inheritance.R > ../R/TL_Toro1_H5_F_carcass_MOI.out

R --vanilla --args ../../data/TL_M_carcass.genes_down_total.txt ../../data/Toro1_M_carcass.genes_down_total.txt ../../data/H6_M_carcass.genes_down_total.txt ../../plots TL_Toro1_H6_M_carcass < ../R/mode_of_inheritance.R > ../R/TL_Toro1_H6_M_carcass_MOI.out

R --vanilla --args ../../data/TL_M_carcass.genes_down_total.txt ../../data/Toro1_M_carcass.genes_down_total.txt ../../data/H5_M_carcass.genes_down_total.txt ../../plots TL_Toro1_H5_M_carcass < ../R/mode_of_inheritance.R > ../R/TL_Toro1_H5_M_carcass_MOI.out

R --vanilla --args ../../data/TL_ovaries.genes_down_total.txt ../../data/Toro1_ovaries.genes_down_total.txt ../../data/H6_ovaries.genes_down_total.txt ../../plots TL_Toro1_H6_ovaries < ../R/mode_of_inheritance.R > ../R/TL_Toro1_H6_ovaries_MOI.out

R --vanilla --args ../../data/TL_ovaries.genes_down_total.txt ../../data/Toro1_ovaries.genes_down_total.txt ../../data/H5_ovaries.genes_down_total.txt ../../plots TL_Toro1_H5_ovaries < ../R/mode_of_inheritance.R > ../R/TL_Toro1_H5_ovaries_MOI.out

R --vanilla --args ../../data/TL_testes.genes_down_total.txt ../../data/Toro1_testes.genes_down_total.txt ../../data/H6_testes.genes_down_total.txt ../../plots TL_Toro1_H6_testes < ../R/mode_of_inheritance.R > ../R/TL_Toro1_H6_testes_MOI.out

R --vanilla --args ../../data/TL_testes.genes_down_total.txt ../../data/Toro1_testes.genes_down_total.txt ../../data/H5_testes.genes_down_total.txt ../../plots TL_Toro1_H5_testes < ../R/mode_of_inheritance.R > ../R/TL_Toro1_H5_testes_MOI.out

mv ../../plots/*.txt ../../data

# combine reg_div and MOI

R --vanilla --args ../../data/TL_Toro1_H6_F_carcass_reg_div.txt ../../data/TL_Toro1_H6_F_carcass_MOI.txt ../../plots TL_Toro1_H6_F_carcass < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H6_F_carcass_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H5_F_carcass_reg_div.txt ../../data/TL_Toro1_H5_F_carcass_MOI.txt ../../plots TL_Toro1_H5_F_carcass < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H5_F_carcass_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H6_M_carcass_reg_div_no_X.txt ../../data/TL_Toro1_H6_M_carcass_MOI.txt ../../plots TL_Toro1_H6_M_carcass < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H6_M_carcass_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H5_M_carcass_reg_div_no_X.txt ../../data/TL_Toro1_H5_M_carcass_MOI.txt ../../plots TL_Toro1_H5_M_carcass < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H5_M_carcass_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H6_ovaries_reg_div.txt ../../data/TL_Toro1_H6_ovaries_MOI.txt ../../plots TL_Toro1_H6_ovaries < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H6_ovaries_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H5_ovaries_reg_div.txt ../../data/TL_Toro1_H5_ovaries_MOI.txt ../../plots TL_Toro1_H5_ovaries < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H5_ovaries_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H6_testes_reg_div_no_X.txt ../../data/TL_Toro1_H6_testes_MOI.txt ../../plots TL_Toro1_H6_testes < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H6_testes_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H5_testes_reg_div_no_X.txt ../../data/TL_Toro1_H5_testes_MOI.txt ../../plots TL_Toro1_H5_testes < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H5_testes_reg_div_MOI.out

# keep the X in male samples

R --vanilla --args ../../data/TL_Toro1_H6_testes.genes_down_ASE.txt ../../plots TL_Toro1_H6_testes < ../R/regulatory_divergence.R > ../R/TL_Toro1_H6_testes_reg_div.out

R --vanilla --args ../../data/TL_Toro1_H5_testes.genes_down_ASE.txt ../../plots TL_Toro1_H5_testes < ../R/regulatory_divergence.R > ../R/TL_Toro1_H5_testes_reg_div.out

R --vanilla --args ../../data/TL_Toro1_H6_M_carcass.genes_down_ASE.txt ../../plots TL_Toro1_H6_M_carcass < ../R/regulatory_divergence.R > ../R/TL_Toro1_H6_M_carcass_reg_div.out

R --vanilla --args ../../data/TL_Toro1_H5_M_carcass.genes_down_ASE.txt ../../plots TL_Toro1_H5_M_carcass < ../R/regulatory_divergence.R > ../R/TL_Toro1_H5_M_carcass_reg_div.out

mv ../../plots/*.txt ../../data
