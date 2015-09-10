data1 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");

data2 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.SNPs.txt",header=TRUE,sep="\t");

data3 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.SNPs.txt",header=TRUE,sep="\t");

data0 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

DGRP <- read.table("/Users/kraigrs/Wittkopp/DGRP/DGRP_line_40_SNPs_const.txt",header=FALSE,sep="\t");

data <- merge(data0,DGRP,by.x=c("chr","pos"),by.y=c("V1","V3"));

nrow(subset(data1,ref_allele/(ref_allele+alt_allele)==0.5))/nrow(data1);
nrow(subset(data2,ref_allele/(ref_allele+alt_allele)==0.5))/nrow(data2);
nrow(subset(data3,ref_allele/(ref_allele+alt_allele)==0.5))/nrow(data3);
nrow(subset(temp,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele)==0.5))/nrow(temp);


### pie charts ###

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_v1_pie.pdf");
pie(c(nrow(subset(data1,ref_allele/(ref_allele+alt_allele)!=0.5)),nrow(subset(data1,ref_allele/(ref_allele+alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();




pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_v2_pie.pdf");
pie(c(nrow(subset(data2,ref_allele/(ref_allele+alt_allele)!=0.5)),nrow(subset(data2,ref_allele/(ref_allele+alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_v3_pie.pdf");
pie(c(nrow(subset(data3,ref_allele/(ref_allele+alt_allele)!=0.5)),nrow(subset(data3,ref_allele/(ref_allele+alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();


pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_v0_pie.pdf");
pie(c(nrow(subset(temp,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele)!=0.5)),nrow(subset(temp,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();

data1 <- subset(data1,ref_allele>0&alt_allele>0);
data2 <- subset(data2,ref_allele>0&alt_allele>0);
data3 <- subset(data3,ref_allele>0&alt_allele>0);

summary(log2(data1$ref_allele/data1$alt_allele));
summary(log2(data2$ref_allele/data2$alt_allele));
summary(log2(data3$ref_allele/data3$alt_allele));

summary(data0$neighbor);
summary(data1$neighbor);
summary(data2$neighbor);
summary(data3$neighbor);


temp1 <- merge(data1,data2,by.x=c("chr","pos"),by.y=c("chr","pos"));
data <- merge(temp1,data3,by.x=c("chr","pos"),by.y=c("chr","pos"));


data <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

ref_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l50_m0.mappability.txt",header=TRUE,sep="\t");
alt_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_alt_line_40.l50_m0.mappability.txt",header=TRUE,sep="\t");
mappability_0mm <- merge(ref_mappability,alt_mappability,by.x=c("chr","position"),by.y=c("chr","position"));

data <- merge(data,mappability_0mm,by.x=c("chr","pos"),by.y=c("chr","position"));

#data <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_alt_line_40.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");

data <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.SNPs.txt",header=TRUE,sep="\t");

###### figure out how many SNPs are within 50b of each SNP ######
### start iterative run here

complete_data <- NULL;

chromosomes <- c("chr2L","chr2R","chr3L","chr3R","chrX");

for(j in 1:length(chromosomes))
{

	arm <- chromosomes[j];

	chrom <- subset(data,chr == arm);
	chrom_sorted <- chrom[order(chrom$pos),];
	chrom <- chrom_sorted;

	steps_right = 0;
	steps_left = 0;

	SNPs_left <- mat.or.vec(nrow(chrom),1);
	SNPs_right <- mat.or.vec(nrow(chrom),1);

	for(i in 1:nrow(chrom))
	{
		if(i == 1)
		{
			SNPs_left[i] <- 0;
			k <- i + 1;
			right <- chrom$pos[k] - chrom$pos[i];
		
			while(right < 50)
			{
				steps_right <- steps_right + 1;
				if(k == nrow(chrom)){break;}
				k <- k + 1;
				right <- chrom$pos[k] - chrom$pos[i];
			}
			SNPs_right[i] = steps_right;
		}
	
		if(i == nrow(chrom))
		{
			SNPs_right[i] <- 0;
			j <- i - 1;
			left <- chrom$pos[i] - chrom$pos[j];
		
			while(left < 50)
			{
				steps_left <- steps_left + 1;
				if(j == 1){break;}
				j <- j - 1;
				left <- chrom$pos[i] - chrom$pos[j];
			}
			SNPs_left[i] = steps_left;
		}
	
		j <- i - 1;
		k <- i + 1;
	
		steps_left <- 0;
		steps_right <- 0;
	
		if(j != 0 & k != (nrow(chrom)+1))
		{
			left <- chrom$pos[i] - chrom$pos[j];
			right <- chrom$pos[k] - chrom$pos[i];
	
			while(left < 50)
			{
				steps_left <- steps_left + 1;
				j <- j - 1;
				if(j == 0){break;}
				left <- chrom$pos[i] - chrom$pos[j];
			}
			SNPs_left[i] = steps_left;
	
			while(right < 50)
			{
				steps_right <- steps_right + 1;
				k <- k + 1;
				if(k == (nrow(chrom)+1)){break;}
				right <- chrom$pos[k] - chrom$pos[i];
			}
		SNPs_right[i] = steps_right;
		}
	}

	neighbor <- SNPs_left+SNPs_right;
	chrom <- cbind(chrom,neighbor);

	complete_data <- rbind(complete_data,chrom);
}
neighbor_plus <- complete_data$neighbor+1;
complete_data <- cbind(complete_data,neighbor_plus);
summary(complete_data$neighbor_plus);

complete_data0 <- complete_data;

# run the previous iteratively to get complete dataset

#write.table(complete_data,file="/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.dm3_ref.bowtie_mm1.neighbor_SNPs.txt",quote=FALSE,sep="\t",row.names=FALSE);

#write.table(ref_bad_SNPs,file="/Users/kraigrs/Wittkopp/Simulations/tiled/ref_bad_SNPs.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE);

#write.table(alt_bad_SNPs,file="/Users/kraigrs/Wittkopp/Simulations/tiled/alt_bad_SNPs.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE);

# 1 mismatch
boxplot(log2(ref_allele.x/alt_allele.x) ~ neighbor,data = subset(complete_data,ref_allele.x > 0 & alt_allele.x > 0),varwidth = TRUE,xlab="Number of neighboring SNPs",ylab="Distribution of log2(ref/alt) among SNPs",border=rgb(0,0,0,0.5),xlim=c(0,13),ylim=c(-6,6));

par(new=TRUE);

# 2 mismatches
boxplot(log2(ref_allele.y/alt_allele.y) ~ neighbor,data = subset(complete_data, ref_allele.y > 0 & alt_allele.y > 0),varwidth = TRUE,xlab="Number of neighboring SNPs",ylab="Distribution of log2(ref/alt)",border=rgb(1,0,0,0.5),xlim=c(0,13),ylim=c(-6,6));

par(new=TRUE);

# 3 mismatches
boxplot(log2(ref_allele/alt_allele) ~ neighbor,data = subset(complete_data, ref_allele > 0 & alt_allele > 0),varwidth = TRUE,xlab="Number of neighboring SNPs",ylab="Distribution of log2(ref/alt)",border=rgb(0,0,1,0.5),xlim=c(0,13),ylim=c(-6,6));


boxplot(log2(ref_allele/alt_allele) ~ neighbor,data = subset(complete_data, ref_allele > 0 & alt_allele > 0),varwidth = TRUE,xlab="Number of neighboring SNPs",ylab="Distribution of log2(ref/alt)",main="DGRP line_40 single reference (3mm)",xlim=c(0,13),ylim=c(-7,7),pars = list(outpch=19,cex=0.2));

par(mfrow=c(2,2));

boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = complete_data, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.2));

abline(h=0.5,lty=2,col="red");

boxplot(dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) ~ neighbor_plus,data = complete_data, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.2));

abline(h=0.5,lty=2,col="red");

#### bubble plot ####

data <- subset(complete_data,is.finite(log2(ref_allele/alt_allele)) & !is.na(log2(ref_allele/alt_allele)));

radii <- sqrt(data$total_overlap/pi);

symbols(
	data$neighbor+1,
	log2(data$ref_allele/data$alt_allele),
	circles=radii);

#### neighbor plus ####

neighbor_plus <- rep(0,nrow(complete_data));

for(i in 1:nrow(complete_data))
{
	if(complete_data$neighbor[i] < 6){neighbor_plus[i] <- complete_data$neighbor[i];}
	else{neighbor_plus[i] <- "6";}
}

temp <- cbind(complete_data,neighbor_plus);

boxplot(log2(ref_allele.x/alt_allele.x) ~ neighbor_plus, data = subset(temp,ref_allele.x > 0 & alt_allele.x > 0),varwidth = TRUE,xlab="Number of neighboring SNPs",ylab="Distribution of log2(ref/alt) among SNPs",main="(# neighboring SNPs capped at 6)",border=rgb(0,0,0,0.5),ylim=c(-3,7));

abline(h = 0,col="red",lty="dashed");

#############################################
# plot medians for the different mismatches #
#############################################

summary(complete_data$neighbor); # range of neighboring SNPs [0,22]

multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

both <- merge(complete_data,multiple,by.x=c("chr","pos"),by.y=c("chr","pos"));

neighbor_plus <- rep(0,nrow(both));
for(i in 1:nrow(both))
{
	if(both$neighbor[i] < 6){neighbor_plus[i] <- both$neighbor[i];}
	else{neighbor_plus[i] <- "6";}
}
both <- cbind(both,neighbor_plus);

#medians <- mat.or.vec(7,5);
#colnames(medians) <- c("neighbor","v1","v2","v3","multiple");

medians <- mat.or.vec(13,4);
colnames(medians) <- c("neighbor","v1","v2","v3");

#vars <- mat.or.vec(7,5);
#colnames(vars) <- c("neighbor","v1","v2","v3","multiple");

vars <- mat.or.vec(13,4);
colnames(vars) <- c("neighbor","v1","v2","v3");

both <- complete_data;
for(i in 1:13)
{
	j <- i-1;
	temp <- subset(both,neighbor == j);
	
	#medians[i,1] <- j;
	#medians[i,2] <- median(log2(temp$ref_allele.x/temp$alt_allele.x));
	#medians[i,3] <- median(log2(temp$ref_allele.y/temp$alt_allele.y));
	#medians[i,4] <- median(log2(temp$ref_allele/temp$alt_allele));
	#medians[i,5] <- median(log2(temp$dm3_ref_ref_allele/temp$dm3_alt_alt_allele));
	
	medians[i,1] <- j;
	medians[i,2] <- median(temp$ref_allele.x/(temp$ref_allele.x+temp$alt_allele.x));
	medians[i,3] <- median(temp$ref_allele.y/(temp$ref_allele.y+temp$alt_allele.y));
	medians[i,4] <- median(temp$ref_allele/(temp$ref_allele+temp$alt_allele));
	#medians[i,5] <- median(temp$dm3_ref_ref_allele/(temp$dm3_ref_ref_allele+temp$dm3_alt_alt_allele));
	
	#vars[i,1] <- j;
	#vars[i,2] <- var(log2(temp$ref_allele.x/temp$alt_allele.x));
	#vars[i,3] <- var(log2(temp$ref_allele.y/temp$alt_allele.y));
	#vars[i,4] <- var(log2(temp$ref_allele/temp$alt_allele));
	#vars[i,5] <- var(log2(temp$dm3_ref_ref_allele/temp$dm3_alt_alt_allele));
	
	vars[i,1] <- j;
	vars[i,2] <- var(temp$ref_allele.x/(temp$ref_allele.x+temp$alt_allele.x));
	vars[i,3] <- var(temp$ref_allele.y/(temp$ref_allele.y+temp$alt_allele.y));
	vars[i,4] <- var(temp$ref_allele/(temp$ref_allele+temp$alt_allele));
	#vars[i,5] <- var(log2(temp$dm3_ref_ref_allele/temp$dm3_alt_alt_allele));
}

plot(medians[,1],medians[,2],type="o",pch=1,col="black",xlim=c(0,12),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="");
par(new=TRUE);
plot(medians[,1],medians[,3],type="o",pch=0,col="black",xlim=c(0,12),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="");
par(new=TRUE);
plot(medians[,1],medians[,4],type="o",pch=2,col="black",xlim=c(0,12),ylim=c(0,1),xlab="# neighboring SNPs",ylab="median ref/(ref+alt)",main="Compare medians");

#legend("bottomright",legend=c("1 mm","2 mm","3 mm"),fill=c("black","red","green"),bty="n");
legend("bottomright",legend=c("1 mm","2 mm","3 mm"),pch=c(1,0,2),bty="n");


plot(vars[,1],vars[,2],type="o",col="black",xlim=c(0,6),ylim=c(0,3),xaxt="n",yaxt="n",xlab="",ylab="");
par(new=TRUE);
plot(vars[,1],vars[,3],type="o",col="red",xlim=c(0,6),ylim=c(0,3),xaxt="n",yaxt="n",xlab="",ylab="");
par(new=TRUE);
plot(vars[,1],vars[,4],type="o",col="green",xlim=c(0,6),ylim=c(0,3),xlab="# neighboring SNPs",ylab="variance log2(ASE)");

legend("topleft",legend=c("1 mismatch","2 mismatches","3 mismatches"),fill=c("black","red","green"));

# compare distributions of log2(ASE) for differing numbers of mismatches

total1 <- subset(complete_data,ref_allele.x > 0 & alt_allele.x > 0);
total2 <- subset(complete_data,ref_allele.y > 0 & alt_allele.y > 0);
total3 <- subset(complete_data,ref_allele > 0 & alt_allele > 0);

obj1 <- hist(log2(total1$ref_allele.x/total1$alt_allele.x),breaks=seq(-4.75,6.5,0.25));
obj2 <- hist(log2(total2$ref_allele.y/total2$alt_allele.y),breaks=seq(-4.75,6.5,0.25));
obj3 <- hist(log2(total3$ref_allele/total3$alt_allele),breaks=seq(-4.75,6.5,0.25));

mat <- cbind(obj1$counts/nrow(total1),obj2$counts/nrow(total2),obj3$counts/nrow(total3));

barplot(t(mat),beside=TRUE,names.arg=seq(-4.50,6.5,0.25),xlab="log2(ref/alt)",ylab="Proportion",main="SNP-based measurements of ASE",col=c("black","red","green"));

legend("topright",legend=c("1 mismatch","2 mismatches","3 mismatches"),fill=c("black","red","green"));


# compare distributions of fraction of reference allele for differing numbers of mismatches

obj1 <- hist(complete_data$ref_allele.x/(complete_data$ref_allele.x+complete_data$alt_allele.x),breaks=20);
obj2 <- hist(complete_data$ref_allele.y/(complete_data$ref_allele.y+complete_data$alt_allele.y),breaks=20);
obj3 <- hist(complete_data$ref_allele/(complete_data$ref_allele+complete_data$alt_allele),breaks=20);

mat <- cbind(obj1$counts/nrow(total1),obj2$counts/nrow(total2),obj3$counts/nrow(total3));

barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="Fraction of reference allele",ylab="Proportion",main="SNP-based measurements of ASE",col=c("black","red","green"));

legend("topright",legend=c("1 mismatch","2 mismatches","3 mismatches"),fill=c("black","red","green"));



temp <- cbind(obj1$breaks[2:length(obj1$breaks)],obj1$counts,
			  obj2$breaks[2:length(obj2$breaks)],obj2$counts,
			  obj3$breaks[2:length(obj3$breaks)],obj3$counts);

temp <- temp[14:nrow(temp),];

par(mfrow=c(1,3));

barplot(temp[,2]/nrow(total1),names.arg=temp[,1],xlab="",ylab="Proportion",main="\n\n1 mismatch",col="gray",ylim=c(0,1),cex.axis=1.5,cex.names=1.5,cex.lab=1.5,cex.main=1.5);
barplot(temp[,4]/nrow(total2),names.arg=temp[,3],xlab="Fraction of reference allele",ylab="",main="SNP-based ASE using single reference\n\n2 mismatches",col="gray",ylim=c(0,1),cex.axis=1.5,cex.names=1.5,cex.lab=1.5,cex.main=1.5);
barplot(temp[,6]/nrow(total3),names.arg=temp[,5],xlab="",ylab="",main="\n\n3 mismatches",col="gray",ylim=c(0,1),cex.axis=1.5,cex.names=1.5,cex.lab=1.5,cex.main=1.5);

pie(c(nrow(subset(total1,log2(ref_allele.x/alt_allele.x) != 0)),nrow(subset(total1,log2(ref_allele.x/alt_allele.x) == 0))),col=c("gray","white"),labels="");

pie(c(nrow(subset(total2,log2(ref_allele.y/alt_allele.y) != 0)),nrow(subset(total2,log2(ref_allele.y/alt_allele.y) == 0))),col=c("gray","white"),labels="");

pie(c(nrow(subset(total3,log2(ref_allele/alt_allele) != 0)),nrow(subset(total3,log2(ref_allele/alt_allele) == 0))),col=c("gray","white"),labels="");

pie(c(nrow(subset(total1,ref_allele.x/(ref_allele.x+alt_allele.x) != 0.5)),nrow(subset(total1,ref_allele.x/(ref_allele.x+alt_allele.x) == 0.5))),col=c("gray","white"),labels=c("AI","no AI"));

pie(c(nrow(subset(total2,ref_allele.y/(ref_allele.y+alt_allele.y) != 0.5)),nrow(subset(total2,ref_allele.y/(ref_allele.y+alt_allele.y) == 0.5))),col=c("gray","white"),labels=c("AI","no AI"));

pie(c(nrow(subset(total3,ref_allele/(ref_allele+alt_allele) != 0.5)),nrow(subset(total3,ref_allele/(ref_allele+alt_allele) == 0.5))),col=c("gray","white"),labels=c("AI","no AI"));

# new plots 

nrow(subset(both,ref_allele.x>0&alt_allele.x>0&ref_allele.y>0&alt_allele.y>0&ref_allele>0&alt_allele>0&dm3_ref_ref_allele>0&dm3_alt_alt_allele>0));

posASE <- subset(both,ref_allele.x>0&alt_allele.x>0&ref_allele.y>0&alt_allele.y>0&ref_allele>0&alt_allele>0&dm3_ref_ref_allele>0&dm3_alt_alt_allele>0);

props <- c(
	nrow( subset(posASE,log2(ref_allele.x/alt_allele.x) != 0) )/nrow(posASE),
	nrow( subset(posASE,log2(ref_allele.y/alt_allele.y) != 0) )/nrow(posASE),
	nrow( subset(posASE,log2(ref_allele/alt_allele) != 0) )/nrow(posASE),
	nrow( subset(posASE,log2(dm3_ref_ref_allele/dm3_alt_alt_allele) != 0) )/nrow(posASE)
);

barplot(props,names.arg=c(1,2,3,0),ylim=c(0,1),xlab="Number of mismatches",ylab="Proportion of AI",main="Comparison of single and multiple genomes");

nrow(subset(both,ref_allele.x>0&alt_allele.x>0&ref_allele.y>0&alt_allele.y>0&ref_allele>0&alt_allele>0&dm3_ref_ref_allele>0&dm3_alt_alt_allele>0));

posASE <- both

props <- c(
	nrow( subset(both,ref_allele.x/(ref_allele.x+alt_allele.x) != 0.5) )/nrow(posASE),
	nrow( subset(both,ref_allele.y/(ref_allele.y+alt_allele.y) != 0.5) )/nrow(posASE),
	nrow( subset(both,ref_allele/(ref_allele+alt_allele) != 0.5) )/nrow(posASE),
	nrow( subset(both,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5) )/nrow(posASE)
);

barplot(props,names.arg=c(1,2,3,0),ylim=c(0,1),xlab="Number of mismatches",ylab="Proportion of AI",main="Comparison of single and multiple genomes");

###################
# compare methods #
###################

##### compare single vs. multiple ASE measurements in exons

single <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.exons.txt",header=TRUE,sep="\t");
multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.exons.txt",header=TRUE,sep="\t");

single_multiple <- merge(single, multiple,by.x="gene_exon",by.y="gene_exon");

posASE <- subset(single_multiple,dm3.x>0&line_40.x>0&dm3.y>0&line_40.y>0);

obj1 <- hist(log2(posASE$dm3.x/posASE$line_40.x),breaks=seq(-6,7,0.25));
obj2 <- hist(log2(posASE$dm3.y/posASE$line_40.y),breaks=seq(-6,7,0.25));
mat <- cbind(obj1$counts/nrow(posASE),obj2$counts/nrow(posASE));

barplot(t(mat),beside=TRUE,names.arg=seq(-4.75,8,0.25),xlab="log2(reference allele/alternative allele)",ylab="Proportion",main="Exon-based measurements of ASE",col=c("black","gray"));

legend("topleft",legend=c("single genome, AI = 0.482","multiple genomes, AI = 0.011"),fill=c("black","gray"),bty="n");

temp <- cbind(obj1$breaks[2:length(obj1$breaks)],obj1$counts,
			  obj2$breaks[2:length(obj2$breaks)],obj2$counts);

#temp <- temp[14:nrow(temp),];
temp <- temp[18:47,];
mat <- cbind(temp[,2]/nrow(posASE),temp[,4]/nrow(posASE));

barplot(t(mat),beside=TRUE,names.arg=temp[,1],xlab="log2(dm3/line_40)",ylab="Proportion",main="Exon-based measurements of ASE",col=c("black","gray"),ylim=c(0,1));

legend("topright",legend=c("single genome, AI = 0.482","multiple genomes, AI = 0.011"),fill=c("black","gray"),bty="n");

# AI

posASE_single <- subset(single_multiple,dm3.x>0&line_40.x>0&log2(dm3.x/line_40.x)!=0);
posASE_multiple <- subset(single_multiple,dm3.y>0&line_40.y>0&log2(dm3.y/line_40.y)!=0);

obj1 <- hist(log2(posASE_single$dm3.x/posASE_single$line_40.x),breaks=seq(-6,7,0.25));
obj2 <- hist(log2(posASE_multiple$dm3.y/posASE_multiple$line_40.y),breaks=seq(-6,7,0.25));
mat <- cbind(obj1$counts/nrow(posASE_single),obj2$counts/nrow(posASE_multiple));

temp <- cbind(obj1$breaks[2:length(obj1$breaks)],obj1$counts,
			  obj2$breaks[2:length(obj2$breaks)],obj2$counts);

#temp <- temp[14:nrow(temp),];
temp <- temp[18:47,];
mat <- cbind(temp[,2]/nrow(posASE_single),temp[,4]/nrow(posASE_multiple));

barplot(t(mat),beside=TRUE,names.arg=temp[,1],xlab="",ylab="",main="",col=c("black","gray"),ylim=c(0,1));

legend("topleft",legend=c("single (n = 31,334)","multiple (n = 674)"),fill=c("black","gray"));
legend("topright",legend="Exons with detectable ASE\nshowing imbalance");

pie(c(nrow(posASE_single),nrow(posASE_multiple)),labels=c("Single","Multiple"),col=c("black","gray"));

##### compare single vs. multiple ASE measurements in SNPs

single <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");
multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

single_multiple <- merge(single,multiple,by.x=c("chr","pos"),by.y=c("chr","pos"));

posASE <- subset(single_multiple,ref_allele>0&alt_allele>0&dm3_ref_ref_allele>0&dm3_alt_alt_allele>0);

obj1 <- hist(log2(posASE$ref_allele/posASE$alt_allele),breaks=seq(-6,7,0.25));
obj2 <- hist(log2(posASE$dm3_ref_ref_allele/posASE$dm3_alt_alt_allele),breaks=seq(-6,7,0.25));
mat <- cbind(obj1$counts/nrow(posASE),obj2$counts/nrow(posASE));

temp <- cbind(obj1$breaks[2:length(obj1$breaks)],obj1$counts,
			  obj2$breaks[2:length(obj2$breaks)],obj2$counts);

#temp <- temp[14:nrow(temp),];
temp <- temp[18:47,];
mat <- cbind(temp[,2]/nrow(posASE),temp[,4]/nrow(posASE));

barplot(t(mat),beside=TRUE,names.arg=temp[,1],xlab="log2(dm3/line_40)",ylab="Proportion",main="SNP-based measurements of ASE",col=c("black","gray"),ylim=c(0,1));

legend("topright",legend=c("single genome, AI = 0.507","multiple genomes, AI = 0.003"),fill=c("black","gray"),bty="n");
legend("topright",legend="78,860 SNPs with detectable ASE");

# AI

posASE_single <- subset(single_multiple,ref_allele>0&alt_allele>0&log2(ref_allele/alt_allele)!=0);
posASE_multiple <- subset(single_multiple,dm3_ref_ref_allele>0&dm3_alt_alt_allele>0&log2(dm3_ref_ref_allele/dm3_alt_alt_allele)!=0);

obj1 <- hist(log2(posASE_single$ref_allele/posASE_single$alt_allele),breaks=seq(-6,7,0.25));
obj2 <- hist(log2(posASE_multiple$dm3_ref_ref_allele/posASE_multiple$dm3_alt_alt_allele),breaks=seq(-6,7,0.25));
mat <- cbind(obj1$counts/nrow(posASE_single),obj2$counts/nrow(posASE_multiple));

temp <- cbind(obj1$breaks[2:length(obj1$breaks)],obj1$counts,
			  obj2$breaks[2:length(obj2$breaks)],obj2$counts);

#temp <- temp[14:nrow(temp),];
temp <- temp[18:47,];
mat <- cbind(temp[,2]/nrow(posASE_single),temp[,4]/nrow(posASE_multiple));

barplot(t(mat),beside=TRUE,names.arg=temp[,1],xlab="",ylab="",main="",col=c("black","gray"),ylim=c(0,1));

legend("topleft",legend=c("single (n = 31,334)","multiple (n = 674)"),fill=c("black","gray"));
legend("topright",legend="Exons with detectable ASE\nshowing imbalance");

pie(c(nrow(posASE_single),nrow(posASE_multiple)),labels=c("Single","Multiple"),col=c("black","gray"));


###### redone using fraction of reference allele
##### compare single vs. multiple ASE measurements in SNPs

single <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");
multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

single_multiple <- merge(single,multiple,by.x=c("chr","pos"),by.y=c("chr","pos"));

posASE <- single_multiple;

obj1 <- hist(posASE$ref_allele/(posASE$ref_allele+posASE$alt_allele),breaks=20);
obj2 <- hist(posASE$dm3_ref_ref_allele/(posASE$dm3_ref_ref_allele+posASE$dm3_alt_alt_allele),breaks=20);
mat <- cbind(obj1$counts/nrow(posASE),obj2$counts/nrow(posASE));

temp <- cbind(obj1$breaks[2:length(obj1$breaks)],obj1$counts,
			  obj2$breaks[2:length(obj2$breaks)],obj2$counts);

#temp <- temp[14:nrow(temp),];
#temp <- temp[18:47,];
mat <- cbind(temp[,2]/nrow(posASE),temp[,4]/nrow(posASE));

barplot(t(mat),beside=TRUE,names.arg=temp[,1],xlab="Fraction of reference allele",ylab="Proportion",main="SNP-based measurements of ASE",col=c("black","gray"),ylim=c(0,1));

legend("topright",legend=c("single genome, AI = 0.507","multiple genomes, AI = 0.003"),fill=c("black","gray"),bty="n");
legend("topright",legend="78,860 SNPs with detectable ASE");

# AI

posASE_single <- subset(single_multiple,ref_allele/(ref_allele+alt_allele)!=0.5);
posASE_multiple <- subset(single_multiple,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele)!=0.5);

obj1 <- hist(posASE_single$ref_allele/(posASE_single$ref_allele+posASE_single$alt_allele),breaks=20);
obj2 <- hist(posASE_multiple$dm3_ref_ref_allele/(posASE_multiple$dm3_ref_ref_allele+posASE_multiple$dm3_alt_alt_allele),breaks=20);
mat <- cbind(obj1$counts/nrow(posASE_single),obj2$counts/nrow(posASE_multiple));

temp <- cbind(obj1$breaks[2:length(obj1$breaks)],obj1$counts,
			  obj2$breaks[2:length(obj2$breaks)],obj2$counts);

#temp <- temp[14:nrow(temp),];
#temp <- temp[18:47,];
mat <- cbind(temp[,2]/nrow(posASE_single),temp[,4]/nrow(posASE_multiple));

barplot(t(mat),beside=TRUE,names.arg=temp[,1],xlab="",ylab="",main="",col=c("black","gray"),ylim=c(0,1));

legend("topleft",legend=c("single (n = 31,334)","multiple (n = 674)"),fill=c("black","gray"));
legend("topright",legend="Exons with detectable ASE\nshowing imbalance");

pie(c(nrow(posASE_single),nrow(posASE_multiple)),labels=c("Single","Multiple"),col=c("black","gray"));

#########################

###############
# mappability #
###############

dm3_ref_exon_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons.dm3_ref.l50_m0.mappability.txt",header=TRUE,sep="\t");
dm3_alt_exon_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons.dm3_alt_line_40.l50_m0.mappability.txt",header=TRUE,sep="\t");

exon_mappability = merge(dm3_ref_exon_mappability,dm3_alt_exon_mappability,by.x="locus",by.y="locus");
dm3_ref_avg <- exon_mappability$sum.x/exon_mappability$length.x;
dm3_alt_avg <- exon_mappability$sum.y/exon_mappability$length.y;

mappability <- cbind(exon_mappability,dm3_ref_avg,dm3_alt_avg);

multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.exons.txt",header=TRUE,sep="\t");

SNPsInExons <- read.table("/Users/kraigrs/Wittkopp/DGRP/DGRP_line_40_SNPs_in_const.txt",header=FALSE,sep="\t");

temp1 <- merge(multiple,mappability,by.x="gene_exon",by.y="locus");
temp2 <- merge(temp1,SNPsInExons,by.x="gene_exon",by.y="V8");

weird <- subset(temp2,dm3>0&line_40>0&V7==0&log2(dm3/line_40)!=0)[,c(1:4,11,12,19)];

AI <- subset(temp2,dm3>0&line_40>0&V7>0&log2(dm3/line_40)!=0);
AI_prop <- nrow(subset(AI,dm3_ref_avg != 1 | dm3_alt_avg != 1))/nrow(AI);
call <- rep("AI",nrow(AI));
AI <- cbind(AI,call);

ASE <- subset(temp2,dm3>0&line_40>0&V7>0&log2(dm3/line_40)==0);
ASE_prop <- nrow(subset(ASE,dm3_ref_avg != 1 | dm3_alt_avg != 1))/nrow(ASE);
call <- rep("ASE",nrow(ASE));
ASE <- cbind(ASE,call);

data <- rbind(AI,ASE);

boxplot(log2(dm3_ref_avg/dm3_alt_avg) ~ call, data=data,ylab="log2(mappability ratio)",varwidth=TRUE);

props <- c(AI_prop,ASE_prop);
barplot(props,names.arg=c("AI\nn = 300","not AI\nn = 23,902"),ylim=c(0,1),ylab="Proportion filtered out in each category",main="Filtering imperfect mappability retains many exons");

barplot(c(ASE_prop,AI_prop),xlim=c(0,1),horiz=TRUE,names.arg=c("No AI\nn = 23,902","AI\nn = 300"),width=0.3,xlab="Proportion removed due to imperfect mappability",main="Exons");

############

# 1 mismatch

mappability_1mm <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l50_m1.mappability.txt",header=TRUE,sep="\t");

data <- merge(complete_data1,mappability_1mm,by.x=c("chr","pos"),by.y=c("chr","position"));
write.table(data,file="/Users/kraigrs/Desktop/exons/sim_SNPs_1mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- subset(data,neighbor<1);
write.table(temp,file="/Users/kraigrs/Desktop/exons/sim_SNPs_1mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

perfect <- subset(temp,sum/length == 1);
write.table(perfect,file="/Users/kraigrs/Desktop/exons/sim_SNPs_1mm_unbiased_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5))*100; # % equal ASRA and imperfect mappability

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length == 1))/nrow(subset(temp,sum/length == 1))*100;

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

nonsig <- subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5);
sig <- subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5);
nrow(subset(sig,sum/length < 1))/nrow(sig)*100;
nrow(subset(nonsig,sum/length < 1))/nrow(nonsig)*100;

nonsig_hist <- hist( (nonsig$sum/nonsig$length) ,breaks=seq(0,1,0.05));
sig_hist <- hist( (sig$sum/sig$length) ,breaks=seq(0,1,0.05));
mat <- cbind(nonsig_hist$counts/nrow(nonsig),sig_hist$counts/nrow(sig));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="",ylab="",main="",col=c("white","grey"));

data <- merge(complete_data1,mappability_1mm,by.x=c("chr","pos"),by.y=c("chr","position"));

perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_1mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,1],ylim=c(0,1),col="black",type="b");
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,2],ylim=c(0,1),col="grey",type="b");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_1mm_pie_by_mapp.pdf");
pie(nrow(perfect),nrow(imperfect),labels="");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_1mm_boxplot_perfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = perfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_1mm_boxplot_imperfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = imperfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

d1_perfect <- density(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),from=0,to=1);
d1_imperfect <- density(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),from=0,to=1);

plot(d1_perfect,xlim=c(0,1),ylim=c(0,10),main="",col="blue",xlab="");
par(new=TRUE);
plot(d1_imperfect,xlim=c(0,1),ylim=c(0,10),main="",col="red",xlab="");

# 2 mismatch

mappability_2mm <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l50_m2.mappability.txt",header=TRUE,sep="\t");

data <- merge(complete_data2,mappability_2mm,by.x=c("chr","pos"),by.y=c("chr","position"));
write.table(data,file="/Users/kraigrs/Desktop/exons/sim_SNPs_2mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- subset(data,neighbor<2);
write.table(temp,file="/Users/kraigrs/Desktop/exons/sim_SNPs_2mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

perfect <- subset(temp,sum/length == 1);
write.table(perfect,file="/Users/kraigrs/Desktop/exons/sim_SNPs_2mm_unbiased_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length == 1))/nrow(subset(temp,sum/length == 1))*100;

nrow(subset(data,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length == 1))/nrow(subset(data,ref_allele/(ref_allele+alt_allele) == 0.5))*100; #  
nrow(subset(data,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length != 1))/nrow(subset(data,ref_allele/(ref_allele+alt_allele) == 0.5))*100; # 

nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5))*100; # % equal ASRA and imperfect mappability

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);


nonsig <- subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5);
sig <- subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5);
nrow(subset(sig,sum/length < 1))/nrow(sig)*100;
nrow(subset(nonsig,sum/length < 1))/nrow(nonsig)*100;

nonsig_hist <- hist( (nonsig$sum/nonsig$length) ,breaks=seq(0,1,0.05));
sig_hist <- hist( (sig$sum/sig$length) ,breaks=seq(0,1,0.05));
mat <- cbind(nonsig_hist$counts/nrow(nonsig),sig_hist$counts/nrow(sig));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="",ylab="",main="",col=c("white","grey"));

data <- merge(complete_data2,mappability_2mm,by.x=c("chr","pos"),by.y=c("chr","position"));

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_2mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,1],ylim=c(0,1),col="black",type="b");
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,2],ylim=c(0,1),col="grey",type="b");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_2mm_boxplot_perfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = perfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_2mm_boxplot_imperfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = imperfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_2mm_pie_by_mapp.pdf");
pie(nrow(perfect),nrow(imperfect),labels="");
dev.off();

d2_perfect <- density(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele));
d2_imperfect <- density(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele));

plot(d2_perfect,xlim=c(0,1),ylim=c(0,10),main="",col="blue",xlab="");
par(new=TRUE);
plot(d2_imperfect,xlim=c(0,1),ylim=c(0,10),main="",col="red",xlab="");

# 3 mismatch

mappability_3mm <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l50_m3.mappability.txt",header=TRUE,sep="\t");

data <- merge(complete_data3,mappability_3mm,by.x=c("chr","pos"),by.y=c("chr","position"));
write.table(data,file="/Users/kraigrs/Desktop/exons/sim_SNPs_3mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- subset(data,neighbor<3);
write.table(temp,file="/Users/kraigrs/Desktop/exons/sim_SNPs_3mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

perfect <- subset(temp,sum/length == 1);
write.table(perfect,file="/Users/kraigrs/Desktop/exons/sim_SNPs_3mm_unbiased_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length == 1))/nrow(subset(temp,sum/length == 1))*100;

nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5))*100; # % equal ASRA and imperfect mappability

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

nonsig <- subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5);
sig <- subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5);
nrow(subset(sig,sum/length < 1))/nrow(sig)*100;
nrow(subset(nonsig,sum/length < 1))/nrow(nonsig)*100;

nonsig_hist <- hist( (nonsig$sum/nonsig$length) ,breaks=seq(0,1,0.05));
sig_hist <- hist( (sig$sum/sig$length) ,breaks=seq(0,1,0.05));
mat <- cbind(nonsig_hist$counts/nrow(nonsig),sig_hist$counts/nrow(sig));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="",ylab="",main="",col=c("white","grey"));

data <- merge(complete_data3,mappability_3mm,by.x=c("chr","pos"),by.y=c("chr","position"));

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_3mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,1],ylim=c(0,1),col="black",type="b");
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,2],ylim=c(0,1),col="grey",type="b");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_3mm_pie_by_mapp.pdf");
pie(nrow(perfect),nrow(imperfect),labels="");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_3mm_boxplot_perfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = perfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_3mm_boxplot_imperfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = imperfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

d3_perfect <- density(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele));
d3_imperfect <- density(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele));

plot(d3_perfect,xlim=c(0,1),ylim=c(0,10),main="",col="blue",xlab="");
par(new=TRUE);
plot(d3_imperfect,xlim=c(0,1),ylim=c(0,10),main="",col="red",xlab="");

# 0 mismatches

ref_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l50_m0.mappability.txt",header=TRUE,sep="\t");
alt_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_alt_line_40.l50_m0.mappability.txt",header=TRUE,sep="\t");
mappability_0mm <- merge(ref_mappability,alt_mappability,by.x=c("chr","position"),by.y=c("chr","position"));

data <- merge(complete_data0,mappability_0mm,by.x=c("chr","pos"),by.y=c("chr","position"));

temp <- data;

write.table(complete_data0,file="/Users/kraigrs/Desktop/exons/sim_SNPs_0mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

write.table(complete_data0,file="/Users/kraigrs/Desktop/exons/sim_SNPs_0mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

perfect <- subset(complete_data0,((sum.x/length.x)+(sum.y/length.y)) == 2);
write.table(perfect,file="/Users/kraigrs/Desktop/exons/sim_SNPs_0mm_unbiased_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) == 0.5 & ((sum.x/length.x)+(sum.y/length.y)) == 2))/nrow(subset(temp,((sum.x/length.x)+(sum.y/length.y)) == 2))*100;


nrow(subset(temp,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5 & ((sum.x/length.x)+(sum.y/length.y)) < 2))/nrow(subset(temp,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) == 0.5 & ((sum.x/length.x)+(sum.y/length.y)) < 2))/nrow(subset(temp,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) == 0.5))*100; # % equal ASRA and imperfect mappability

perfect <- subset(temp,((sum.x/length.x)+(sum.y/length.y)) == 2);
imperfect <- subset(temp,((sum.x/length.x)+(sum.y/length.y)) < 2);

nonsig <- subset(temp,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) == 0.5);
sig <- subset(temp,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5);
#nrow(subset(sig,sum.x/length.x + sum.y/length.y != 2))/nrow(sig)*100;
#nrow(subset(nonsig,sum.x/length.x + sum.y/length.y != 2))/nrow(nonsig)*100;

nonsig_hist <- hist( (nonsig$sum.x/nonsig$length.x)+(nonsig$sum.y/nonsig$length.y) ,breaks=seq(0,2,0.1));
sig_hist <- hist( (sig$sum.x/sig$length.x)+(sig$sum.y/sig$length.y) ,breaks=seq(0,2,0.1));
mat <- cbind(nonsig_hist$counts/nrow(nonsig),sig_hist$counts/nrow(sig));
barplot(t(mat),beside=TRUE,names.arg=seq(0.1,2,0.1),xlab="",ylab="",main="",col=c("white","grey"));

plot(log2((data$sum.x/data$length.x)/(data$sum.y/data$length.y)),data$dm3_ref_ref_allele/(data$dm3_ref_ref_allele+data$dm3_alt_alt_allele),xlab="",ylab="",pch=19,col=rgb(0,0,0,0.3),cex=0.7);


perfect <- subset(temp,((sum.x/length.x)+(sum.y/length.y)) == 2);
imperfect <- subset(temp,((sum.x/length.x)+(sum.y/length.y)) != 2);

perfect_hist <- hist(perfect$dm3_ref_ref_allele/(perfect$dm3_ref_ref_allele+perfect$dm3_alt_alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$dm3_ref_ref_allele/(imperfect$dm3_ref_ref_allele+imperfect$dm3_alt_alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_0mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,1],ylim=c(0,1),col="black",type="b");
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,2],ylim=c(0,1),col="grey",type="b");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_0mm_pie_by_mapp.pdf");
pie(nrow(perfect),nrow(imperfect),labels="");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_0mm_boxplot_perfect.pdf");
boxplot(dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) ~ neighbor_plus,data = perfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_0mm_boxplot_imperfect.pdf");
boxplot(dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) ~ neighbor_plus,data = imperfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

d0_perfect <- density(perfect$dm3_ref_ref_allele/(perfect$dm3_ref_ref_allele+perfect$dm3_alt_alt_allele));
d0_imperfect <- density(imperfect$dm3_ref_ref_allele/(imperfect$dm3_ref_ref_allele+imperfect$dm3_alt_alt_allele));

plot(d0_perfect,xlim=c(0,1),ylim=c(0,10),main="",col="blue",xlab="");
par(new=TRUE);
plot(d0_imperfect,xlim=c(0,1),ylim=c(0,10),main="",col="red",xlab="");


#########################

plot((data$sum.x/data$length.x)/(data$sum.x/data$length.x+data$sum.y/data$length.y),data$dm3_ref_ref_allele/(data$dm3_ref_ref_allele+data$dm3_alt_alt_allele),xlab="mappability",ylab="proportion of reference allele",pch=19,col=rgb(0,0,0,0.5),cex=0.5)

data <- subset(data,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5 & (sum.x/length.x)/(sum.x/length.x+sum.y/length.y) != 0.5);

data <- subset(data,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5);

nrow(subset(data, sign((sum.x/length.x)/(sum.x/length.x+sum.y/length.y)-0.5) == sign(dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele)-0.5) ));

multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

data <- merge(multiple,SNP_mappability,by.x=c("chr","pos"),by.y=c("chr","position"));

nrow(subset(data,dm3_ref_ref_allele > 0 & dm3_alt_alt_allele > 0 & log2(dm3_ref_ref_allele/dm3_alt_alt_allele) != 0));
nrow(subset(data,dm3_ref_ref_allele > 0 & dm3_alt_alt_allele > 0 & log2(dm3_ref_ref_allele/dm3_alt_alt_allele) != 0 & log2((sum.x/length.x)/(sum.y/length.y)) != 0));

nrow(subset(data,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5));
nrow(subset(data,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5 & log2((sum.x/length.x)/(sum.y/length.y)) != 0));

poor_map <- subset(data,dm3_ref_ref_allele > 0 & dm3_alt_alt_allele > 0 & log2(dm3_ref_ref_allele/dm3_alt_alt_allele) != 0 & log2((sum.x/length.x)/(sum.y/length.y)) != 0);
nrow(subset(poor_map,sign(log2(dm3_ref_ref_allele/dm3_alt_alt_allele)) == sign(log2((sum.x/length.x)/(sum.y/length.y)))));

AI <- subset(data,ref_allele.x > 0 & alt_allele.x > 0 & log2(ref_allele.x/alt_allele.x) != 0);
AI_prop <- nrow(subset(AI,sum/length != 1))/nrow(AI);

ASE <- subset(data,ref_allele.x > 0 & alt_allele.x > 0 & log2(ref_allele.x/alt_allele.x) == 0);
ASE_prop <- nrow(subset(ASE,sum/length != 1))/nrow(ASE);

# SNPs with no neighbors and AI

temp1 <- subset(AI,neighbor==0);
AI_prop <- nrow(subset(temp1,sum/length != 1))/nrow(temp1); 

temp2 <- subset(ASE,neighbor==0);
ASE_prop <- nrow(subset(temp2,sum/length != 1))/nrow(temp2); 

props <- c(AI_prop,ASE_prop);
barplot(props,names.arg=c("AI\nn = 179","not AI\nn = 38,819"),ylim=c(0,1),ylab="Proportion filtered out in each category",main="Filtering imperfect mappability retains many SNPs");

barplot(c(ASE_prop,AI_prop),xlim=c(0,1),horiz=TRUE,names.arg=c("No AI\nn = 179","AI\nn = 38,819"),width=0.3,xlab="Proportion removed due to imperfect mappability",main="SNPs");

##################################
# plot mappability across genome #
##################################

library(lattice);

dm3_ref_exon_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons.dm3_ref.l50_m0.mappability.txt",header=TRUE,sep="\t");

dm3_alt_exon_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons.dm3_alt_line_40.l50_m0.mappability.txt",header=TRUE,sep="\t");

exon_mappability = merge(dm3_ref_exon_mappability,dm3_alt_exon_mappability,by.x="locus",by.y="locus");
dm3_ref_avg <- exon_mappability$sum.x/exon_mappability$length.x;
dm3_alt_avg <- exon_mappability$sum.y/exon_mappability$length.y;

mappability <- cbind(exon_mappability,dm3_ref_avg,dm3_alt_avg);

multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.exons.txt",header=TRUE,sep="\t");

SNPsInExons <- read.table("/Users/kraigrs/Wittkopp/DGRP/DGRP_line_40_SNPs_in_const.txt",header=FALSE,sep="\t");

temp1 <- merge(multiple,mappability,by.x="gene_exon",by.y="locus");
temp2 <- merge(temp1,SNPsInExons,by.x="gene_exon",by.y="V8");

temp3 <- subset(temp2,dm3 > 0 & line_40 > 0);

xyplot(dm3_ref_avg/dm3_alt_avg ~ (V2+V3)/2 | V1, data = temp3,xlab="Midpoint of exon",ylab="dm3/line_40 mappability",main="Mappability across both allele-specific genomes",pch=19,cex=0.4,col=rgb(0,0,0,0.3),layout=c(2,3));

xyplot(dm3_ref_avg/dm3_alt_avg ~ log2(dm3/line_40) | V1, data = temp3,xlab="log2(dm3/line_40) ASE",ylab="dm3/line_40 mappability",main="Concordance of mappability direction and bias in ASE",pch=19,cex=0.4,col=rgb(0,0,0,0.3));

# how many exons show opposite signs of mappability and bias? 0.2%, so 99.8% show same sign, awesome!
# also, of the 660 that do show differential mappability, 237 favor dm3_ref and 423 favor dm3_alt
sum( sign(log2(temp3$dm3_ref_avg/temp3$dm3_alt_avg)) != sign(log2(temp3$dm3/temp3$line_40)) );

# correlation? r^2 = 0.6
cor(log2(temp3$dm3_ref_avg/temp3$dm3_alt_avg),log2(temp3$dm3_ref/temp3$dm3_alt))^2;

###############

SNPsInExons <- read.table("/Users/kraigrs/Wittkopp/DGRP/DGRP_line_40_SNPs_in_const.txt",header=FALSE,sep="\t");

v0 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.exons.txt",header=TRUE,sep="\t");
v1 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.exons.txt",header=TRUE,sep="\t");
v2 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.exons.txt",header=TRUE,sep="\t");
v3 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.exons.txt",header=TRUE,sep="\t");

v0 <- merge(v0,SNPsInExons,by.x="gene_exon",by.y="V8");
v1 <- merge(v1,SNPsInExons,by.x="gene_exon",by.y="V8");
v2 <- merge(v2,SNPsInExons,by.x="gene_exon",by.y="V8");
v3 <- merge(v3,SNPsInExons,by.x="gene_exon",by.y="V8");

nrow(subset(v0,dm3>0&line_40>0&log2(dm3/line_40)!=0&V7>0))/nrow(v0);
nrow(subset(v1,dm3>0&line_40>0&log2(dm3/line_40)!=0&V7>0))/nrow(v1);
nrow(subset(v2,dm3>0&line_40>0&log2(dm3/line_40)!=0&V7>0))/nrow(v2);
nrow(subset(v3,dm3>0&line_40>0&log2(dm3/line_40)!=0&V7>0))/nrow(v3);

nrow(subset(v0,dm3>0&line_40>0&V7>0));
nrow(subset(v1,dm3>0&line_40>0&V7>0));
nrow(subset(v2,dm3>0&line_40>0&V7>0));
nrow(subset(v3,dm3>0&line_40>0&V7>0));


# mappability for single-genome approach

data <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");

SNP_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l50_m1.mappability.txt",header=TRUE,sep="\t");

data <- merge(complete_data,SNP_mappability,by.x=c("chr","pos"),by.y=c("chr","position"));

AI <- subset(data, ref_allele/(ref_allele+alt_allele) != 0.5));
AI_prop <- nrow(subset(AI,sum/length != 1))/nrow(AI);

ASE <- subset(data,ref_allele.x > 0 & alt_allele.x > 0 & log2(ref_allele.x/alt_allele.x) == 0);
ASE_prop <- nrow(subset(ASE,sum/length != 1))/nrow(ASE);

# SNPs with no neighbors and AI

temp1 <- subset(AI,neighbor==0);
AI_prop <- nrow(subset(temp1,sum/length != 1))/nrow(temp1); 

temp2 <- subset(ASE,neighbor==0);
ASE_prop <- nrow(subset(temp2,sum/length != 1))/nrow(temp2);

##############
# alternative allele used as genome

data <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_alt_line_40.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");
neighbors_plus <- data$neighbors+1;
data <- cbind(data,neighbors_plus);

boxplot(alt_allele/(alt_allele+ref_allele) ~ neighbors_plus,data = data, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.2));

abline(h=0.5,lty=2,col="red");

pie(c(nrow(subset(data,alt_allele/(ref_allele+alt_allele)!=0.5)),nrow(data)),col=c("gray","white"),labels="");
