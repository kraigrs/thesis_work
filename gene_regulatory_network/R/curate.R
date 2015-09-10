zhrXz30 <- read.table("~/Wittkopp/gene_regulatory_network/data/zhrXz30_mRNA_meta_gene_output_121312_classified_C.txt",header=TRUE);

data <- zhrXz30[,c(1:3,5,13,14,16,36:50)]

colnames(data) <- c("gene","P.s1","P.s2","P.total","H.s1","H.s2","H.total","cis","trans","allcis","alltrans","conserved1","compensatory","ambiguous","cisplustrans","cisXtrans","additive","underdominant","overdominant","P1dominant","P2dominant","conserved2")

write.table(data,file="~/Wittkopp/gene_regulatory_network/data/zhrXz30_regulatory_divergence_classification.txt",quote=FALSE,sep="\t",row.names=FALSE)

z30Xzhr <- read.table("~/Wittkopp/gene_regulatory_network/data/z30Xzhr_mRNA_meta_gene_output_121312_classified_C.txt",header=TRUE);

data <- z30Xzhr[,c(1:3,5,13,14,16,36:50)]

colnames(data) <- c("gene","P.s1","P.s2","P.total","H.s1","H.s2","H.total","cis","trans","allcis","alltrans","conserved1","compensatory","ambiguous","cisplustrans","cisXtrans","additive","underdominant","overdominant","P1dominant","P2dominant","conserved2")

write.table(data,file="~/Wittkopp/gene_regulatory_network/data/z30Xzhr_regulatory_divergence_classification.txt",quote=FALSE,sep="\t",row.names=FALSE)