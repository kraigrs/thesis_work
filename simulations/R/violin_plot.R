library(vioplot);

data1 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");

relevant <- read.table("/Users/kraigrs/Wittkopp/Simulations/sim_regs.36b_windows_SNPs.neighbors.txt",header=TRUE,sep="\t");

data1 <- merge(data1,relevant,by.x=c("chr","pos"),by.y=c("chr","pos"));

summary(data1);

temp1 <- subset(data1,relevant == 1);
temp2 <- subset(data1,relevant == 2);
temp3 <- subset(data1,relevant == 3);
temp4 <- subset(data1,relevant == 4);
temp5 <- subset(data1,relevant == 5);
temp6 <- subset(data1,relevant == 6);
temp7 <- subset(data1,relevant == 7);
temp8 <- subset(data1,relevant == 8);
temp9 <- subset(data1,relevant == 9);

v1 <- temp1$ref_allele/(temp1$ref_allele+temp1$alt_allele);
v2 <- temp2$ref_allele/(temp2$ref_allele+temp2$alt_allele);
v3 <- temp3$ref_allele/(temp3$ref_allele+temp3$alt_allele);
v4 <- temp4$ref_allele/(temp4$ref_allele+temp4$alt_allele);
v5 <- temp5$ref_allele/(temp5$ref_allele+temp5$alt_allele);
v6 <- temp6$ref_allele/(temp6$ref_allele+temp6$alt_allele);
v7 <- temp7$ref_allele/(temp7$ref_allele+temp7$alt_allele);
v8 <- temp8$ref_allele/(temp8$ref_allele+temp8$alt_allele);
v9 <- temp9$ref_allele/(temp9$ref_allele+temp9$alt_allele);

v1 <- v1[!is.na(v1)];
v2 <- v2[!is.na(v2)];
v3 <- v3[!is.na(v3)];
v4 <- v4[!is.na(v4)];
v5 <- v5[!is.na(v5)];
v6 <- v6[!is.na(v6)];
v7 <- v7[!is.na(v7)];
v8 <- v8[!is.na(v8)];
v9 <- v9[!is.na(v9)];

vioplot(v1, v2, v3, v4, v5, v6, v7, v8, v9,
		names = c(1:9),
		col="cyan");
abline(h=0.5,lty=2,col="red");