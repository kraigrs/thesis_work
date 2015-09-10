#!/bin/sh

R --vanilla --args ../../data/TL_Toro1_H6_F_carcass_reg_div.txt ../../data/TL_Toro1_H6_F_carcass_MOI.txt ../../plots TL_Toro1_H6_F_carcass < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H6_F_carcass_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H5_F_carcass_reg_div.txt ../../data/TL_Toro1_H5_F_carcass_MOI.txt ../../plots TL_Toro1_H5_F_carcass < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H5_F_carcass_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H6_M_carcass_reg_div_no_X.txt ../../data/TL_Toro1_H6_M_carcass_MOI.txt ../../plots TL_Toro1_H6_M_carcass < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H6_M_carcass_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H5_M_carcass_reg_div_no_X.txt ../../data/TL_Toro1_H5_M_carcass_MOI.txt ../../plots TL_Toro1_H5_M_carcass < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H5_M_carcass_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H6_ovaries_reg_div.txt ../../data/TL_Toro1_H6_ovaries_MOI.txt ../../plots TL_Toro1_H6_ovaries < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H6_ovaries_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H5_ovaries_reg_div.txt ../../data/TL_Toro1_H5_ovaries_MOI.txt ../../plots TL_Toro1_H5_ovaries < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H5_ovaries_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H6_testes_reg_div_no_X.txt ../../data/TL_Toro1_H6_testes_MOI.txt ../../plots TL_Toro1_H6_testes < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H6_testes_reg_div_MOI.out

R --vanilla --args ../../data/TL_Toro1_H5_testes_reg_div_no_X.txt ../../data/TL_Toro1_H5_testes_MOI.txt ../../plots TL_Toro1_H5_testes < ../R/regulatory_divergence_MOI.R > ../R/TL_Toro1_H5_testes_reg_div_MOI.out