library(venneuler);

# compare SAMtools and GATK

# single nucleotide variants

# Dpse_TL
Dpse_TL_SNPs <- venneuler(c(SAMtools=19679,GATK=775905,"SAMtools&GATK"=635307));

# Dbog_Toro1
Dbog_Toro1_SNPs <- venneuler(c(SAMtools=52979,GATK=920311,"SAMtools&GATK"=993016));

par(mfrow=c(1,2));
plot(Dpse_TL_SNPs);
plot(Dbog_Toro1_SNPs);




library(VennDiagram);

venn.diagram(list(SAMtools=1:19679,GATK=))