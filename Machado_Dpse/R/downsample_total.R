library(BiasedUrn)

set.seed(12345);

data1 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_F_carcass.genes.txt",header=TRUE,sep="\t");
data2 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_M_carcass.genes.txt",header=TRUE,sep="\t");
data3 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_ovaries.genes.txt",header=TRUE,sep="\t");
data4 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_testes.genes.txt",header=TRUE,sep="\t");

data5 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_F_carcass.genes.txt",header=TRUE,sep="\t");
data6 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_M_carcass.genes.txt",header=TRUE,sep="\t");
data7 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_ovaries.genes.txt",header=TRUE,sep="\t");
data8 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_testes.genes.txt",header=TRUE,sep="\t");

data9 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_F_carcass.genes.txt",header=TRUE,sep="\t");
data10 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_M_carcass.genes.txt",header=TRUE,sep="\t");
data11 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_ovaries.genes.txt",header=TRUE,sep="\t");
data12 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_testes.genes.txt",header=TRUE,sep="\t");

data13 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_F_carcass.genes.txt",header=TRUE,sep="\t");
data14 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_M_carcass.genes.txt",header=TRUE,sep="\t");
data15 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_ovaries.genes.txt",header=TRUE,sep="\t");
data16 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_testes.genes.txt",header=TRUE,sep="\t");

vec1 <- as.vector(t(data1[,c(2:4)]));
vec2 <- as.vector(t(data2[,c(2:4)]));
vec3 <- as.vector(t(data3[,c(2:4)]));
vec4 <- as.vector(t(data4[,c(2:4)]));
vec5 <- as.vector(t(data5[,c(2:4)]));
vec6 <- as.vector(t(data6[,c(2:4)]));
vec7 <- as.vector(t(data7[,c(2:4)]));
vec8 <- as.vector(t(data8[,c(2:4)]));
vec9 <- as.vector(t(data9[,c(2:4)]));
vec10 <- as.vector(t(data10[,c(2:4)]));
vec11 <- as.vector(t(data11[,c(2:4)]));
vec12 <- as.vector(t(data12[,c(2:4)]));
vec13 <- as.vector(t(data13[,c(2:4)]));
vec14 <- as.vector(t(data14[,c(2:4)]));
vec15 <- as.vector(t(data15[,c(2:4)]));
vec16 <- as.vector(t(data16[,c(2:4)]));

mapping <- c(sum(vec1),sum(vec2),sum(vec3),sum(vec4),sum(vec5),sum(vec6),sum(vec7),sum(vec8),
             sum(vec9),sum(vec10),sum(vec11),sum(vec12),sum(vec13),sum(vec14),sum(vec15),sum(vec16));

min_mapping <- min(mapping);

vec1_down <- rMFNCHypergeo(1, vec1, min_mapping, odds=rep(1,length(vec1)), precision = 1E-7);
vec2_down <- rMFNCHypergeo(1, vec2, min_mapping, odds=rep(1,length(vec2)), precision = 1E-7);
vec3_down <- rMFNCHypergeo(1, vec3, min_mapping, odds=rep(1,length(vec3)), precision = 1E-7);
vec4_down <- rMFNCHypergeo(1, vec4, min_mapping, odds=rep(1,length(vec4)), precision = 1E-7);
vec5_down <- rMFNCHypergeo(1, vec5, min_mapping, odds=rep(1,length(vec5)), precision = 1E-7);
vec6_down <- rMFNCHypergeo(1, vec6, min_mapping, odds=rep(1,length(vec6)), precision = 1E-7);
vec7_down <- rMFNCHypergeo(1, vec7, min_mapping, odds=rep(1,length(vec7)), precision = 1E-7);
vec8_down <- rMFNCHypergeo(1, vec8, min_mapping, odds=rep(1,length(vec8)), precision = 1E-7);
vec9_down <- rMFNCHypergeo(1, vec9, min_mapping, odds=rep(1,length(vec9)), precision = 1E-7);
vec10_down <- rMFNCHypergeo(1, vec10, min_mapping, odds=rep(1,length(vec10)), precision = 1E-7);
vec11_down <- rMFNCHypergeo(1, vec11, min_mapping, odds=rep(1,length(vec11)), precision = 1E-7);
vec12_down <- rMFNCHypergeo(1, vec12, min_mapping, odds=rep(1,length(vec12)), precision = 1E-7);
vec13_down <- rMFNCHypergeo(1, vec13, min_mapping, odds=rep(1,length(vec13)), precision = 1E-7);
vec14_down <- rMFNCHypergeo(1, vec14, min_mapping, odds=rep(1,length(vec14)), precision = 1E-7);
vec15_down <- rMFNCHypergeo(1, vec15, min_mapping, odds=rep(1,length(vec15)), precision = 1E-7);
vec16_down <- rMFNCHypergeo(1, vec16, min_mapping, odds=rep(1,length(vec16)), precision = 1E-7);

mat1 <- matrix(vec1_down,nrow=nrow(data1),ncol=3,byrow=TRUE);
mat2 <- matrix(vec2_down,nrow=nrow(data2),ncol=3,byrow=TRUE);
mat3 <- matrix(vec3_down,nrow=nrow(data3),ncol=3,byrow=TRUE);
mat4 <- matrix(vec4_down,nrow=nrow(data4),ncol=3,byrow=TRUE);
mat5 <- matrix(vec5_down,nrow=nrow(data5),ncol=3,byrow=TRUE);
mat6 <- matrix(vec6_down,nrow=nrow(data6),ncol=3,byrow=TRUE);
mat7 <- matrix(vec7_down,nrow=nrow(data7),ncol=3,byrow=TRUE);
mat8 <- matrix(vec8_down,nrow=nrow(data8),ncol=3,byrow=TRUE);
mat9 <- matrix(vec9_down,nrow=nrow(data9),ncol=3,byrow=TRUE);
mat10 <- matrix(vec10_down,nrow=nrow(data10),ncol=3,byrow=TRUE);
mat11 <- matrix(vec11_down,nrow=nrow(data11),ncol=3,byrow=TRUE);
mat12 <- matrix(vec12_down,nrow=nrow(data12),ncol=3,byrow=TRUE);
mat13 <- matrix(vec13_down,nrow=nrow(data13),ncol=3,byrow=TRUE);
mat14 <- matrix(vec14_down,nrow=nrow(data14),ncol=3,byrow=TRUE);
mat15 <- matrix(vec15_down,nrow=nrow(data15),ncol=3,byrow=TRUE);
mat16 <- matrix(vec16_down,nrow=nrow(data16),ncol=3,byrow=TRUE);

genes1 <- data1[,1]; data1_out <- data.frame(genes1,mat1); colnames(data1_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes2 <- data2[,1]; data2_out <- data.frame(genes2,mat2); colnames(data2_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes3 <- data3[,1]; data3_out <- data.frame(genes3,mat3); colnames(data3_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes4 <- data4[,1]; data4_out <- data.frame(genes4,mat4); colnames(data4_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes5 <- data5[,1]; data5_out <- data.frame(genes5,mat5); colnames(data5_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes6 <- data6[,1]; data6_out <- data.frame(genes6,mat6); colnames(data6_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes7 <- data7[,1]; data7_out <- data.frame(genes7,mat7); colnames(data7_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes8 <- data8[,1]; data8_out <- data.frame(genes8,mat8); colnames(data8_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes9 <- data9[,1]; data9_out <- data.frame(genes9,mat9); colnames(data9_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes10 <- data10[,1]; data10_out <- data.frame(genes10,mat10); colnames(data10_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes11 <- data11[,1]; data11_out <- data.frame(genes11,mat11); colnames(data11_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes12 <- data12[,1]; data12_out <- data.frame(genes12,mat12); colnames(data12_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes13 <- data13[,1]; data13_out <- data.frame(genes13,mat13); colnames(data13_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes14 <- data14[,1]; data14_out <- data.frame(genes14,mat14); colnames(data14_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes15 <- data15[,1]; data15_out <- data.frame(genes15,mat15); colnames(data15_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");
genes16 <- data16[,1]; data16_out <- data.frame(genes16,mat16); colnames(data16_out) <- c("gene","Dpse_TL","Dbog_Toro1","Both");

write.table(data1_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_F_carcass.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(data2_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_M_carcass.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(data3_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_ovaries.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(data4_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_testes.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);

write.table(data5_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_F_carcass.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(data6_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_M_carcass.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(data7_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_ovaries.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(data8_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_testes.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);

write.table(data9_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_F_carcass.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(data10_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_M_carcass.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(data11_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_ovaries.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(data12_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_testes.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);

write.table(data13_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_F_carcass.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(data14_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_M_carcass.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(data15_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_ovaries.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(data16_out,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_testes.genes_down_total.txt",quote=FALSE,sep="\t",row.names=FALSE);
