exons <- read.table("/Users/kraigrs/Wittkopp/Simulations/SNPs_in_const.txt",header=FALSE,sep="\t");
lengths <- exons[,3]-exons[,2];
exons <- cbind(exons,lengths);

exprn <- read.table("/Users/kraigrs/Wittkopp/Simulations/zhr_z30_exons_expression.txt",sep="\t");

exons <- merge(exons,exprn,by.x="V8",by.y="V1");

tiled_exons <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.bowtie.exons.txt",header=TRUE,sep="\t");

equal_allele_exons <- read.table("/Users/kraigrs/Wittkopp/Simulations/equal_allele/constExons_single_bp50_error0_equal_allele.bowtie.exons.txt",header=TRUE,sep="\t");

equal_total_exons <- read.table("/Users/kraigrs/Wittkopp/Simulations/equal_total/constExons_single_bp50_error0_equal_total.bowtie.exons.txt",header=TRUE,sep="\t");

#############################
# choose a group to look at #
#############################

merged <- merge(tiled_exons,exons,by.x="gene_exon",by.y="V8");

ASE <- subset(merged,merged$dm3_ref > 0 & merged$dm3_alt > 0 & log2(merged$dm3_ref/merged$dm3_alt) == 0);
no_ASE <- subset(merged,merged$dm3_ref == 0 & merged$dm3_alt == 0 & merged$Both > 0);

AI <- subset(merged,merged$dm3_ref > 0 & merged$dm3_alt > 0 & log2(merged$dm3_ref/merged$dm3_alt) != 0);

expressed <- subset(merged,merged$V2.y > 0);

##############
# make plots #
##############

# proportion of TUs with ASE
nrow(ASE)/nrow(merged);
nASE <- nrow(ASE);

# proportion of TUs without ASE but with measurable expression
nrow(no_ASE)/nrow(merged);

# proportion of TUs displaying allelic imbalance
nrow(AI)/nrow(merged);
nAI <- nrow(AI);

# barplot comparing number of SNPs between TU with and without allelic imbalance

obj1 <- hist(ASE$V7/(ASE$lengths/1000),breaks=seq(0,180,5)); # 240 is based on the maximum between each of ASE and AI, which is 236
obj2 <- hist(AI$V7/(AI$lengths/1000),breaks=seq(0,180,5));
mat <- cbind(obj1$counts/nrow(ASE),obj2$counts/nrow(AI));

obj3 <- hist(log2(AI$dm3_ref/AI$dm3_alt),breaks=40);

par( mfrow = c( 1, 2 ) );

barplot(t(mat),beside=TRUE,names.arg=seq(0,175,5),xlab="Number of SNPs/kb",ylab="Proportion of exons",xlim=c(0,60),main="",col=c("black","gray"));

plot(obj3$breaks[1:length(obj3$breaks)-1],obj3$counts/nrow(AI),type="h",col="black",ylim=c(0,0.3),xlim=c(-7,5),xlab="log2(ref/alt)",ylab="Proportion of exons",main="");


#################

exons <- read.table("/Users/kraigrs/Wittkopp/Simulations/SNPs_in_const.txt",header=FALSE,sep="\t");
lengths <- exons[,3]-exons[,2];
exons <- cbind(exons,lengths);

equal_allele_exons <- read.table("/Users/kraigrs/Wittkopp/Simulations/equal_allele/constExons_single_bp50_error0_equal_allele.bowtie.exons.txt",header=TRUE,sep="\t");
merged1 <- merge(equal_allele_exons,exons,by.x="gene_exon",by.y="V8");

equal_total_exons <- read.table("/Users/kraigrs/Wittkopp/Simulations/equal_total/constExons_single_bp50_error0_equal_total.bowtie.exons.txt",header=TRUE,sep="\t");
merged2 <- merge(equal_total_exons,exons,by.x="gene_exon",by.y="V8");

exprn <- read.table("/Users/kraigrs/Wittkopp/Simulations/zhr_z30_exons_expression.txt",sep="\t");

data1 <- merge(merged1,exprn,by.x="gene_exon",by.y="V1");
data2 <- merge(merged2,exprn,by.x="gene_exon",by.y="V1");

par(mfrow=c(1,2));
plot(log2(data1$dm3_ref/data1$dm3_alt),log2(data1$V2.y),pch=19,col=rgb(0,0,0,0.2),cex=0.3,
	main="Equal allele approach",ylab="log2(number of generated reads)",xlab="log2(ref/alt)");
abline(v = 0,col="red");

plot(log2(data2$dm3_ref/data2$dm3_alt),log2(data2$dm3_ref+data2$dm3_alt+data2$Both),pch=19,col=rgb(0,0,0,0.2),cex=0.3,
	main="Equal total approach",ylab="log2(number of generated reads)",xlab="log2(ref/alt)");
abline(v = 0,col="red");

########################
# binomial exact tests #
########################

cut = 0.05;

tmp1 <- subset(tiled_exons,tiled_exons$dm3_ref > 0 & tiled_exons$dm3_alt > 0);
pvals <- NULL;
for(i in 1:nrow(tmp1))
{
	test <- binom.test(tmp1$dm3_ref[i],tmp1$dm3_ref[i]+tmp1$dm3_alt[i], p = 0.5, alternative = "two.sided", conf.level = 0.95);
	pvals <- c(pvals,test$p.value);
}

FPR <- sum(pvals<cut)/length(pvals);

tmp2 <- subset(equal_allele_exons,equal_allele_exons$dm3_ref > 0 & equal_allele_exons$dm3_alt > 0);
pvals <- NULL;
for(i in 1:nrow(tmp2))
{
	test <- binom.test(tmp2$dm3_ref[i],tmp2$dm3_ref[i]+tmp2$dm3_alt[i], p = 0.5, alternative = "two.sided", conf.level = 0.95);
	pvals <- c(pvals,test$p.value);
}

FPR <- sum(pvals<cut)/length(pvals);

tmp3 <- subset(equal_total_exons,equal_total_exons$dm3_ref > 0 & equal_total_exons$dm3_alt > 0);
pvals <- NULL;
for(i in 1:nrow(tmp3))
{
	test <- binom.test(tmp3$dm3_ref[i],tmp3$dm3_ref[i]+tmp3$dm3_alt[i], p = 0.5, alternative = "two.sided", conf.level = 0.95);
	pvals <- c(pvals,test$p.value);
}
sum(pvals<cut);
length(pvals);
sum(pvals<cut)/length(pvals);

##############
# pie charts #
##############

pie(c(nrow(ASE),nrow(AI)),labels=c("ASE","AI"),col=c("black","gray"));



