#########################
### Create RPKM table ###
#########################

library(edgeR)

counts = read.csv("row_expression_counts_Cavalli.tsv",sep="\t")
rownames(counts)=counts$X
counts$X=NULL

y=DGEList(counts=counts[,c(1,2)], genes=data.frame(length=counts[,3]))
y = calcNormFactors(y)
RPKM = as.data.frame(rpkm(y, log=T))

RPKM$Average = apply(RPKM, 1, mean)

write.table(RPKM,
            'RPKM_expression_AllGenes_mRNA_mESC_Cavalli.tsv',
            sep="\t",
            quote=F)

