# choose a line from DGRP

data <- read.table("/Users/kraigrs/Wittkopp/DGRP/lines_all_chrs.txt",header=TRUE,sep="\t");

# max. homozygous genotypes
data[which.max(data$A+data$C+data$G+data$T),];

# min. heterozygous calls
data[which.min(data$R+data$Y+data$M+data$K+data$S+data$W),];
