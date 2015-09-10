library(BiasedUrn)

set.seed(67891);

data1 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_F_carcass.genes_down_total.txt",header=TRUE,sep="\t");
data2 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_F_carcass.genes_down_total.txt",header=TRUE,sep="\t");

data3 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_M_carcass.genes_down_total.txt",header=TRUE,sep="\t");
data4 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_M_carcass.genes_down_total.txt",header=TRUE,sep="\t");

data5 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_ovaries.genes_down_total.txt",header=TRUE,sep="\t");
data6 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_ovaries.genes_down_total.txt",header=TRUE,sep="\t");

data7 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_testes.genes_down_total.txt",header=TRUE,sep="\t");
data8 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_testes.genes_down_total.txt",header=TRUE,sep="\t");

data9 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_F_carcass.genes_down_total.txt",header=TRUE,sep="\t");
data10 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_M_carcass.genes_down_total.txt",header=TRUE,sep="\t");
data11 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_ovaries.genes_down_total.txt",header=TRUE,sep="\t");
data12 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_testes.genes_down_total.txt",header=TRUE,sep="\t");

data13 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_F_carcass.genes_down_total.txt",header=TRUE,sep="\t");
data14 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_M_carcass.genes_down_total.txt",header=TRUE,sep="\t");
data15 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_ovaries.genes_down_total.txt",header=TRUE,sep="\t");
data16 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_testes.genes_down_total.txt",header=TRUE,sep="\t");

genes <- data1[,1];

genes_12 <- data1[,1];
data_12 <- data.frame(genes_12,data1[,2]+data2[,2],data1[,3]+data2[,3]);
colnames(data_12) <- c("gene","Dpse_TL","Dbog_Toro1");

genes34 <- data3[,1];
data34 <- data.frame(genes34,data3[,2]+data4[,2],data3[,3]+data4[,3]);
colnames(data34) <- c("gene","Dpse_TL","Dbog_Toro1");

genes56 <- data5[,1];
data56 <- data.frame(genes56,data5[,2]+data6[,2],data5[,3]+data6[,3]);
colnames(data56) <- c("gene","Dpse_TL","Dbog_Toro1");

genes78 <- data7[,1];
data78 <- data.frame(genes78,data7[,2]+data8[,2],data7[,3]+data8[,3]);
colnames(data78) <- c("gene","Dpse_TL","Dbog_Toro1");
             
mat1 <- cbind(data_12[,2]+data_12[,3],data9[,2]+data9[,3]); # TL_Toro1_H5_F_carcass
mat2 <- cbind(data_12[,2]+data_12[,3],data13[,2]+data13[,3]); # TL_Toro1_H6_F_carcass

mat3 <- cbind(data34[,2]+data34[,3],data10[,2]+data10[,3]); # TL_Toro1_H5_M_carcass
mat4 <- cbind(data34[,2]+data34[,3],data14[,2]+data14[,3]); # TL_Toro1_H6_M_carcass

mat5 <- cbind(data56[,2]+data56[,3],data11[,2]+data11[,3]); # TL_Toro1_H5_ovaries
mat6 <- cbind(data56[,2]+data56[,3],data15[,2]+data15[,3]); # TL_Toro1_H6_ovaries

mat7 <- cbind(data78[,2]+data78[,3],data12[,2]+data12[,3]); # TL_Toro1_H5_testes
mat8 <- cbind(data78[,2]+data78[,3],data16[,2]+data16[,3]); # TL_Toro1_H6_testes

mat9 <- cbind(data13[,2]+data13[,3],data9[,2]+data9[,3]); # H6_H5_F_carcass
mat10 <- cbind(data14[,2]+data14[,3],data10[,2]+data10[,3]); # H6_H5_M_carcass
mat11 <- cbind(data15[,2]+data15[,3],data11[,2]+data11[,3]); # H6_H5_ovaries
mat12 <- cbind(data16[,2]+data16[,3],data12[,2]+data12[,3]); # H6_H5_F_testes

minima1 <- apply(mat1,1,min);
minima2 <- apply(mat2,1,min);
minima3 <- apply(mat3,1,min);
minima4 <- apply(mat4,1,min);
minima5 <- apply(mat5,1,min);
minima6 <- apply(mat6,1,min);
minima7 <- apply(mat7,1,min);
minima8 <- apply(mat8,1,min);
minima9 <- apply(mat9,1,min);
minima10 <- apply(mat10,1,min);
minima11 <- apply(mat11,1,min);
minima12 <- apply(mat12,1,min);

downsample <- function(x)
{
	tmp1 <- rhyper(1,x[1],x[2],x[5]);
	tmp2 <- rhyper(1,x[3],x[4],x[5]);
	value1 <- c(tmp1,x[5]-tmp1);
	value2 <- c(tmp2,x[5]-tmp2);
	values <- c(value1,value2);
	return(values);
}

temp1 <- t(apply(as.matrix(cbind(data_12[,c(2,3)],data9[,c(2,3)],minima1)),1,downsample));
down1 <- data.frame(genes,temp1); colnames(down1) <- c("gene","Dpse_TL","Dbog_Toro1","Hyb_Dpse","Hyb_Dbog");

temp2 <- t(apply(as.matrix(cbind(data_12[,c(2,3)],data13[,c(2,3)],minima2)),1,downsample));
down2 <- data.frame(genes,temp2); colnames(down2) <- c("gene","Dpse_TL","Dbog_Toro1","Hyb_Dpse","Hyb_Dbog");

temp3 <- t(apply(as.matrix(cbind(data34[,c(2,3)],data10[,c(2,3)],minima3)),1,downsample));
down3 <- data.frame(genes,temp3); colnames(down3) <- c("gene","Dpse_TL","Dbog_Toro1","Hyb_Dpse","Hyb_Dbog");

temp4 <- t(apply(as.matrix(cbind(data34[,c(2,3)],data14[,c(2,3)],minima4)),1,downsample));
down4 <- data.frame(genes,temp4); colnames(down4) <- c("gene","Dpse_TL","Dbog_Toro1","Hyb_Dpse","Hyb_Dbog");

temp5 <- t(apply(as.matrix(cbind(data56[,c(2,3)],data11[,c(2,3)],minima5)),1,downsample));
down5 <- data.frame(genes,temp5); colnames(down5) <- c("gene","Dpse_TL","Dbog_Toro1","Hyb_Dpse","Hyb_Dbog");

temp6 <- t(apply(as.matrix(cbind(data56[,c(2,3)],data15[,c(2,3)],minima6)),1,downsample));
down6 <- data.frame(genes,temp6); colnames(down6) <- c("gene","Dpse_TL","Dbog_Toro1","Hyb_Dpse","Hyb_Dbog");

temp7 <- t(apply(as.matrix(cbind(data78[,c(2,3)],data12[,c(2,3)],minima7)),1,downsample));
down7 <- data.frame(genes,temp7); colnames(down7) <- c("gene","Dpse_TL","Dbog_Toro1","Hyb_Dpse","Hyb_Dbog");

temp8 <- t(apply(as.matrix(cbind(data78[,c(2,3)],data16[,c(2,3)],minima8)),1,downsample));
down8 <- data.frame(genes,temp8); colnames(down8) <- c("gene","Dpse_TL","Dbog_Toro1","Hyb_Dpse","Hyb_Dbog");

temp9 <- t(apply(as.matrix(cbind(data13[,c(2,3)],data9[,c(2,3)],minima9)),1,downsample));
down9 <- data.frame(genes,temp9); colnames(down9) <- c("gene","H6_Dpse","H6_Dbog","H5_Dpse","H5_Dbog");

temp10 <- t(apply(as.matrix(cbind(data14[,c(2,3)],data10[,c(2,3)],minima10)),1,downsample));
down10 <- data.frame(genes,temp10); colnames(down10) <- c("gene","H6_Dpse","H6_Dbog","H5_Dpse","H5_Dbog");

temp11 <- t(apply(as.matrix(cbind(data15[,c(2,3)],data11[,c(2,3)],minima11)),1,downsample));
down11 <- data.frame(genes,temp11); colnames(down11) <- c("gene","H6_Dpse","H6_Dbog","H5_Dpse","H5_Dbog");

temp12 <- t(apply(as.matrix(cbind(data16[,c(2,3)],data12[,c(2,3)],minima12)),1,downsample));
down12 <- data.frame(genes,temp12); colnames(down12) <- c("gene","H6_Dpse","H6_Dbog","H5_Dpse","H5_Dbog");

write.table(down1,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H5_F_carcass.genes_down_ASE.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(down2,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass.genes_down_ASE.txt",quote=FALSE,sep="\t",row.names=FALSE);

write.table(down3,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H5_M_carcass.genes_down_ASE.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(down4,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_M_carcass.genes_down_ASE.txt",quote=FALSE,sep="\t",row.names=FALSE);

write.table(down5,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H5_ovaries.genes_down_ASE.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(down6,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_ovaries.genes_down_ASE.txt",quote=FALSE,sep="\t",row.names=FALSE);

write.table(down7,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H5_testes.genes_down_ASE.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(down8,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_testes.genes_down_ASE.txt",quote=FALSE,sep="\t",row.names=FALSE);

write.table(down9,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_H5_F_carcass.genes_down_ASE.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(down10,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_H5_M_carcass.genes_down_ASE.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(down11,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_H5_ovaries.genes_down_ASE.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(down12,"/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_H5_testes.genes_down_ASE.txt",quote=FALSE,sep="\t",row.names=FALSE);
