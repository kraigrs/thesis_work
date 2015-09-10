############################################################################################
# script to explore different properties of the data
############################################################################################

annotation <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.genes.bed",header=FALSE,sep="\t");
colnames(annotation) <- c("chromosome","start","stop","gene","empty","strand");

# TL_F_carcass

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_F_carcass.genes.txt",header=TRUE,sep="\t");
data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");
data <- data_annotation;

prop <- round(sum(data$Dpse_TL)/(sum(data$Dpse_TL)+sum(data$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_F_carcass.genes.maternal_allele.pdf");

hist(data$Dpse_TL/(data$Dpse_TL+data$Dbog_Toro1),breaks=20,xlab="proportion of maternal allele",main="TL_F_carcass.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");

dev.off();

# TL_ovaries

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_ovaries.genes.txt",header=TRUE,sep="\t");
data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");
data <- data_annotation;

prop <- round(sum(data$Dpse_TL)/(sum(data$Dpse_TL)+sum(data$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_ovaries.genes.maternal_allele.pdf");

hist(data$Dpse_TL/(data$Dpse_TL+data$Dbog_Toro1),breaks=20,xlab="proportion of maternal allele",main="TL_ovaries.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");

dev.off();

# TL_M_carcass

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_M_carcass.genes.txt",header=TRUE,sep="\t");
data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");
data <- data_annotation;

prop <- round(sum(data$Dpse_TL)/(sum(data$Dpse_TL)+sum(data$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_M_carcass.genes.maternal_allele.pdf");

hist(data$Dpse_TL/(data$Dpse_TL+data$Dbog_Toro1),breaks=20,xlab="proportion of maternal allele",main="TL_M_carcass.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");

dev.off();

# TL_testes

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_testes.genes.txt",header=TRUE,sep="\t");
data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");
data <- data_annotation;

prop <- round(sum(data$Dpse_TL)/(sum(data$Dpse_TL)+sum(data$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_testes.genes.maternal_allele.pdf");

hist(data$Dpse_TL/(data$Dpse_TL+data$Dbog_Toro1),breaks=20,xlab="proportion of maternal allele",main="TL_testes.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");

dev.off();

# Toro1_F_carcass

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_F_carcass.genes.txt",header=TRUE,sep="\t");
data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");
data <- data_annotation;

prop <- round(sum(data$Dbog_Toro1)/(sum(data$Dpse_TL)+sum(data$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/Toro1_F_carcass.genes.maternal_allele.pdf");

hist(data$Dbog_Toro1/(data$Dpse_TL+data$Dbog_Toro1),breaks=20,xlab="proportion of maternal allele",main="Toro1_F_carcass.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");

dev.off();

# Toro1_ovaries

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_ovaries.genes.txt",header=TRUE,sep="\t");
data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");
data <- data_annotation;

prop <- round(sum(data$Dbog_Toro1)/(sum(data$Dpse_TL)+sum(data$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/Toro1_ovaries.genes.maternal_allele.pdf");

hist(data$Dbog_Toro1/(data$Dpse_TL+data$Dbog_Toro1),breaks=20,xlab="proportion of maternal allele",main="Toro1_ovaries.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");

dev.off();

# Toro1_M_carcass

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_M_carcass.genes.txt",header=TRUE,sep="\t");
data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");
data <- data_annotation;

prop <- round(sum(data$Dbog_Toro1)/(sum(data$Dpse_TL)+sum(data$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/Toro1_M_carcass.genes.maternal_allele.pdf");

hist(data$Dbog_Toro1/(data$Dpse_TL+data$Dbog_Toro1),breaks=20,xlab="proportion of maternal allele",main="Toro1_M_carcass.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");

dev.off();

# Toro1_testes

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_testes.genes.txt",header=TRUE,sep="\t");
data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");
data <- data_annotation;

prop <- round(sum(data$Dbog_Toro1)/(sum(data$Dpse_TL)+sum(data$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/Toro1_testes.genes.maternal_allele.pdf");

hist(data$Dbog_Toro1/(data$Dpse_TL+data$Dbog_Toro1),breaks=20,xlab="proportion of maternal allele",main="Toro1_testes.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");

dev.off();

# H6_M_carcass

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_M_carcass.genes.txt",header=TRUE,sep="\t");
data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");
data_X <- data_annotation[grep("^X",data_annotation$chromosome,perl=TRUE),];

prop <- round(sum(data_X$Dpse_TL)/(sum(data_X$Dpse_TL)+sum(data_X$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/H6_M_carcass.genes.maternal_allele.pdf");

hist(data_X$Dpse_TL/(data_X$Dpse_TL+data_X$Dbog_Toro1),breaks=20,xlab="proportion of maternal allele",main="H6_M_carcass.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");

dev.off();

# H6_testes

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_testes.genes.txt",header=TRUE,sep="\t");
data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");
data_X <- data_annotation[grep("^X",data_annotation$chromosome,perl=TRUE),];

prop <- round(sum(data_X$Dpse_TL)/(sum(data_X$Dpse_TL)+sum(data_X$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/H6_testes.genes.maternal_allele.pdf");

hist(data_X$Dpse_TL/(data_X$Dpse_TL+data_X$Dbog_Toro1),breaks=20,xlab="proportion of maternal allele",main="H6_testes.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");

dev.off();

# H5_M_carcass

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_M_carcass.genes.txt",header=TRUE,sep="\t");
data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");
data_X <- data_annotation[grep("^X",data_annotation$chromosome,perl=TRUE),];

prop <- round(sum(data_X$Dbog_Toro1)/(sum(data_X$Dpse_TL)+sum(data_X$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/H5_M_carcass.genes.maternal_allele.pdf");

hist(data_X$Dbog_Toro1/(data_X$Dpse_TL+data_X$Dbog_Toro1),breaks=20,xlab="proportion of maternal allele",main="H5_M_carcass.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");

dev.off();

# H5_testes

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H5_testes.genes.txt",header=TRUE,sep="\t");
data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");
data_X <- data_annotation[grep("^X",data_annotation$chromosome,perl=TRUE),];

prop <- round(sum(data_X$Dbog_Toro1)/(sum(data_X$Dpse_TL)+sum(data_X$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/H5_testes.genes.maternal_allele.pdf");

hist(data_X$Dbog_Toro1/(data_X$Dpse_TL+data_X$Dbog_Toro1),breaks=20,xlab="proportion of maternal allele",main="H5_testes.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");

dev.off();
