setwd("/Users/kraigrs/Wittkopp/RegEvol_head_data/genomes/coverage");

zhr <- read.table("zhr.coverage.txt",sep="\t");

hist(zhr$V6/zhr$V7);