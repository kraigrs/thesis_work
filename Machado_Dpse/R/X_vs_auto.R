############################################################################################
# script to compare XL (X) and XR (3L)
############################################################################################

library(vioplot);

annotation <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.genes.bed",header=FALSE,sep="\t");
colnames(annotation) <- c("chromosome","start","stop","gene","empty","strand");

# F_carcass

reg_div <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass_reg_div.txt",header=TRUE,sep="\t");
MOI <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass_MOI.txt",header=TRUE,sep="\t");
data <- merge(reg_div,MOI,by.x=1,by.y=1);

data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");

data_XL <- data_annotation[grep("^XL",data_annotation$chromosome,perl=TRUE),];
data_XR <- data_annotation[grep("^XR",data_annotation$chromosome,perl=TRUE),];
data_auto <- data_annotation[grep("^[^X]",data_annotation$chromosome,perl=TRUE),];

temp <- c(sum(is.finite(data_XL$percent_cis)),sum(is.finite(data_XR$percent_cis)),sum(is.finite(data_auto$percent_cis)));
parms <- scale(sqrt(temp),center=FALSE);

plot(0:4,rep(1,5),type='n',xlim=c(0,4),ylim=c(0,1),ylab="",xlab="",xaxt="n",lty=3,main="%cis");

if(sum(is.finite(data_XL$percent_cis)) > 0){vioplot(data_XL$percent_cis[is.finite(data_XL$percent_cis)],
	at=1,wex=parms[1],add=TRUE,names="XL",col="gray",pchMed=19);}

if(sum(is.finite(data_XR$percent_cis)) > 0){vioplot(data_XR$percent_cis[is.finite(data_XR$percent_cis)],
	at=2,wex=parms[2],add=TRUE,names="XR",col="gray",pchMed=19);}

if(sum(is.finite(data_auto$percent_cis)) > 0){vioplot(data_auto$percent_cis[is.finite(data_auto$percent_cis)],
	at=3,wex=parms[3],add=TRUE,names="autosomes",col="gray",pchMed=19);}
	
reg_div <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H5_F_carcass_reg_div.txt",header=TRUE,sep="\t");
MOI <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H5_F_carcass_MOI.txt",header=TRUE,sep="\t");
data <- merge(reg_div,MOI,by.x=1,by.y=1);

data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");

data_XL <- data_annotation[grep("^XL",data_annotation$chromosome,perl=TRUE),];
data_XR <- data_annotation[grep("^XR",data_annotation$chromosome,perl=TRUE),];
data_auto <- data_annotation[grep("^[^X]",data_annotation$chromosome,perl=TRUE),];

temp <- c(sum(is.finite(data_XL$percent_cis)),sum(is.finite(data_XR$percent_cis)),sum(is.finite(data_auto$percent_cis)));
parms <- scale(sqrt(temp),center=FALSE);

plot(0:4,rep(1,5),type='n',xlim=c(0,4),ylim=c(0,1),ylab="",xlab="",xaxt="n",lty=3,main="%cis");

if(sum(is.finite(data_XL$percent_cis)) > 0){vioplot(data_XL$percent_cis[is.finite(data_XL$percent_cis)],
	at=1,wex=parms[1],add=TRUE,names="XL",col="gray",pchMed=19);}

if(sum(is.finite(data_XR$percent_cis)) > 0){vioplot(data_XR$percent_cis[is.finite(data_XR$percent_cis)],
	at=2,wex=parms[2],add=TRUE,names="XR",col="gray",pchMed=19);}

if(sum(is.finite(data_auto$percent_cis)) > 0){vioplot(data_auto$percent_cis[is.finite(data_auto$percent_cis)],
	at=3,wex=parms[3],add=TRUE,names="autosomes",col="gray",pchMed=19);}
	
# ovaries

reg_div <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_ovaries_reg_div.txt",header=TRUE,sep="\t");
MOI <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_ovaries_MOI.txt",header=TRUE,sep="\t");
data <- merge(reg_div,MOI,by.x=1,by.y=1);

data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");

data_XL <- data_annotation[grep("^XL",data_annotation$chromosome,perl=TRUE),];
data_XR <- data_annotation[grep("^XR",data_annotation$chromosome,perl=TRUE),];
data_auto <- data_annotation[grep("^[^X]",data_annotation$chromosome,perl=TRUE),];

temp <- c(sum(is.finite(data_XL$percent_cis)),sum(is.finite(data_XR$percent_cis)),sum(is.finite(data_auto$percent_cis)));
parms <- scale(sqrt(temp),center=FALSE);

plot(0:4,rep(1,5),type='n',xlim=c(0,4),ylim=c(0,1),ylab="",xlab="",xaxt="n",lty=3,main="%cis");

if(sum(is.finite(data_XL$percent_cis)) > 0){vioplot(data_XL$percent_cis[is.finite(data_XL$percent_cis)],
	at=1,wex=parms[1],add=TRUE,names="XL",col="gray",pchMed=19);}

if(sum(is.finite(data_XR$percent_cis)) > 0){vioplot(data_XR$percent_cis[is.finite(data_XR$percent_cis)],
	at=2,wex=parms[2],add=TRUE,names="XR",col="gray",pchMed=19);}

if(sum(is.finite(data_auto$percent_cis)) > 0){vioplot(data_auto$percent_cis[is.finite(data_auto$percent_cis)],
	at=3,wex=parms[3],add=TRUE,names="autosomes",col="gray",pchMed=19);}
	
cor(data_XL$par1,data_XL$par2,use="pairwise.complete.obs",method="spearman");
cor(data_XR$par1,data_XR$par2,use="pairwise.complete.obs",method="spearman");
cor(data_auto$par1,data_auto$par2,use="pairwise.complete.obs",method="spearman");

t.test(data_XL$percent_cis[is.finite(data_XL$percent_cis)],data_XR$percent_cis[is.finite(data_XR$percent_cis)],alternative="two.sided");
	
reg_div <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H5_ovaries_reg_div.txt",header=TRUE,sep="\t");
MOI <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H5_ovaries_MOI.txt",header=TRUE,sep="\t");
data <- merge(reg_div,MOI,by.x=1,by.y=1);

data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");

data_XL <- data_annotation[grep("^XL",data_annotation$chromosome,perl=TRUE),];
data_XR <- data_annotation[grep("^XR",data_annotation$chromosome,perl=TRUE),];
data_auto <- data_annotation[grep("^[^X]",data_annotation$chromosome,perl=TRUE),];

temp <- c(sum(is.finite(data_XL$percent_cis)),sum(is.finite(data_XR$percent_cis)),sum(is.finite(data_auto$percent_cis)));
parms <- scale(sqrt(temp),center=FALSE);

plot(0:4,rep(1,5),type='n',xlim=c(0,4),ylim=c(0,1),ylab="",xlab="",xaxt="n",lty=3,main="%cis");

if(sum(is.finite(data_XL$percent_cis)) > 0){vioplot(data_XL$percent_cis[is.finite(data_XL$percent_cis)],
	at=1,wex=parms[1],add=TRUE,names="XL",col="gray",pchMed=19);}

if(sum(is.finite(data_XR$percent_cis)) > 0){vioplot(data_XR$percent_cis[is.finite(data_XR$percent_cis)],
	at=2,wex=parms[2],add=TRUE,names="XR",col="gray",pchMed=19);}

if(sum(is.finite(data_auto$percent_cis)) > 0){vioplot(data_auto$percent_cis[is.finite(data_auto$percent_cis)],
	at=3,wex=parms[3],add=TRUE,names="autosomes",col="gray",pchMed=19);}
	
	
	
