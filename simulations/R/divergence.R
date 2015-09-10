genome <- 23011544 + 21146708 + 24543557 + 27905053 + 22422827;
SNPs <- 602328;

SNPsperKb <- SNPs/(genome/1000);

constExons <- read.table("/Users/kraigrs/Wittkopp/Simulations/sim_regs.bed",sep="\t");

SNPsInExons <- 376355;

coding_seq <- sum(constExons$V3-constExons$V2);

