############################################################################################
# script to characterize mode of inheritance
############################################################################################

#Usage: R --vanilla --args parent1 parent2 hybrid < mode_of_inheritance.R > mode_of_inheritance.out

############################################################################################
# functions

# mode of inheritance
mode_of_inheritance <- function(locus)
{
	p1 <- locus[1];
	p2 <- locus[2];
	h <- locus[3];
	
	if((p1/h < 0.8 | p1/h > 1.25) & (p2/h >= 0.8 & p2/h <= 1.25)){return("par2dominant");}
	else if((p1/h >= 0.8 & p1/h <= 1.25) & (p2/h < 0.8 | p2/h > 1.25)){return("par1dominant");}
	else if((p1/h < 0.8 | p1/h > 1.25) & (p2/h < 0.8 | p2/h > 1.25) & h > min(p1,p2) & h < max(p1,p2)){return("additive");}
	else if((p1/h < 0.8 | p1/h > 1.25) & (p2/h < 0.8 | p2/h > 1.25) & h < min(p1,p2)){return("underdominant");}
	else if((p1/h < 0.8 | p1/h > 1.25) & (p2/h < 0.8 | p2/h > 1.25) & h > max(p1,p2)){return("overdominant");}
	else if((p1/h >= 0.8 & p1/h <= 1.25) & (p2/h >= 0.8 & p2/h <= 1.25)){return("similar");}
}

############################################################################################
# read in the data, apply cutoff(s) and filters, perform functions

args <- commandArgs(trailingOnly = TRUE);
#print(args);
parent1 <- read.table(args[1],header=TRUE,sep="\t");
parent2 <- read.table(args[2],header=TRUE,sep="\t");
hybrid <- read.table(args[3],header=TRUE,sep="\t");
dir <- args[4];
out <- args[5];

#parent1 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_F_carcass.genes_down_total.txt",header=TRUE,sep="\t");
#parent2 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Toro1_F_carcass.genes_down_total.txt",header=TRUE,sep="\t");
#hybrid <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/H6_F_carcass.genes_down_total.txt",header=TRUE,sep="\t");
#dir <- "/Users/kraigrs/Wittkopp/Machado_Dpse/plots";
#out <- "TL_Toro1_H6_carcass.genes";

cutoff <- 20;

reliable <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/reliable_genes.txt",header=TRUE,sep="\t");

temp <- merge(parent1,parent2,by.x=1,by.y=1);
data <- merge(temp,hybrid,by.x=1,by.y=1);

temp <- merge(data,reliable,by.x="gene",by.y="gene");
data <- temp;

data_cutoff <- data[which(data[2]+data[3]+data[4] >= cutoff & 
                          data[5]+data[6]+data[7] >= cutoff &
                          data[8]+data[9]+data[10] >= cutoff),];

MOI <- apply(as.matrix(cbind(apply(data_cutoff[,c(2:4)],1,sum),
							 apply(data_cutoff[,c(5:7)],1,sum),
							 apply(data_cutoff[,c(8:10)],1,sum))),1,mode_of_inheritance);

data_final <- cbind(data_cutoff,MOI);

matrix <- as.matrix(cbind(apply(data_cutoff[,c(2:4)],1,sum),
				apply(data_cutoff[,c(5:7)],1,sum),
				apply(data_cutoff[,c(8:10)],1,sum)));
				
gene <- data_final[,1];
table <- data.frame(gene,matrix,MOI);
colnames(table) <- c("gene","par1","par2","hyb","MOI");
write.table(table,file=paste(dir,"/",out,"_MOI.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE);
							 
#axis <- ceiling(max(c(max(abs(log2(matrix[,3])-log2(matrix[,1]))),max(abs(log2(matrix[,3])-log2(matrix[,2]))))));
axis <- 5;

similar <- data_final[data_final$MOI=="similar",];
par1dominant <- data_final[data_final$MOI=="par1dominant",];
par2dominant <- data_final[data_final$MOI=="par2dominant",];
additive <- data_final[data_final$MOI=="additive",];
underdominant <- data_final[data_final$MOI=="underdominant",];
overdominant <- data_final[data_final$MOI=="overdominant",];

pdf(file=paste(dir,"/",out,"_MOI.pdf",sep=""));

plot(log2(similar[,8]+similar[,9]+similar[,10])-log2(similar[,2]+similar[,3]+similar[,4]),
     log2(similar[,8]+similar[,9]+similar[,10])-log2(similar[,5]+similar[,6]+similar[,7]),
     xlim=c(-axis,axis),ylim=c(-axis,axis),col="cyan",pch=19,cex=0.25,
     xlab="log2(hybrid)-log2(Dpse)",ylab="log2(hybrid)-log2(Dbog)",main=paste(out,sep=""));

points(log2(par1dominant[,8]+par1dominant[,9]+par1dominant[,10])-log2(par1dominant[,2]+par1dominant[,3]+par1dominant[,4]),
       log2(par1dominant[,8]+par1dominant[,9]+par1dominant[,10])-log2(par1dominant[,5]+par1dominant[,6]+par1dominant[,7]),
       xlim=c(-axis,axis),ylim=c(-axis,axis),col="green",pch=19,cex=0.25);

points(log2(par2dominant[,8]+par2dominant[,9]+par2dominant[,10])-log2(par2dominant[,2]+par2dominant[,3]+par2dominant[,4]),
       log2(par2dominant[,8]+par2dominant[,9]+par2dominant[,10])-log2(par2dominant[,5]+par2dominant[,6]+par2dominant[,7]),
       xlim=c(-axis,axis),ylim=c(-axis,axis),col="magenta",pch=19,cex=0.25);
       
points(log2(additive[,8]+additive[,9]+additive[,10])-log2(additive[,2]+additive[,3]+additive[,4]),
       log2(additive[,8]+additive[,9]+additive[,10])-log2(additive[,5]+additive[,6]+additive[,7]),
       xlim=c(-axis,axis),ylim=c(-axis,axis),col="orange",pch=19,cex=0.25);
       
points(log2(underdominant[,8]+underdominant[,9]+underdominant[,10])-log2(underdominant[,2]+underdominant[,3]+underdominant[,4]),
       log2(underdominant[,8]+underdominant[,9]+underdominant[,10])-log2(underdominant[,5]+underdominant[,6]+underdominant[,7]),
       xlim=c(-axis,axis),ylim=c(-axis,axis),col="red",pch=19,cex=0.25);

points(log2(overdominant[,8]+overdominant[,9]+overdominant[,10])-log2(overdominant[,2]+overdominant[,3]+overdominant[,4]),
       log2(overdominant[,8]+overdominant[,9]+overdominant[,10])-log2(overdominant[,5]+overdominant[,6]+overdominant[,7]),
       xlim=c(-axis,axis),ylim=c(-axis,axis),col="blue",pch=19,cex=0.25);

dev.off();


data_final$MOI <- factor(data_final$MOI,levels=c("similar","par1dominant","par2dominant","additive","underdominant","overdominant"));

pdf(file=paste(dir,"/",out,"_MOI_pie.pdf",sep=""));
pie(table(data_final$MOI),labels="",col=c("cyan","green","magenta","orange","red","blue"));
dev.off();
