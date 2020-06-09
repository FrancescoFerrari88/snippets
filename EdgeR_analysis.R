# EdgeR analysis

library(edgeR)

setwd("/data/manke/group/ferrari/PhD_project/reference_datasets/Ferrari_iNPC_DMSOvsEPZ_RNA-Seq/downstream_analysis/EdgeR/")

### LOAD COUNT TABLE ###
count=read.csv("counts.tsv", sep='\t', header = T, comment.char = '#')
rownames(count) = count$X
count$X = NULL
colnames(count) = c("NPC48h_DMSO_2","NPC48h_EPZ_2a",
                    "NPC48h_DMSO_4","NPC48h_EPZ_4a",
                    "NPC48h_DMSO_5","NPC48h_EPZ_5a")
count = count[,c(1,3,5,2,4,6)]

#######################
### mRNA processing ###
#######################


### prepare design table ###
design_table=data.frame(row.names = colnames(count),
                        condition = c(rep("DMSO",3),rep("EPZ",3)),
                        mESC_passage = c("P5","P4", "P7","P5","P4","P7"))

all(rownames(design_table) == colnames(count))

y = DGEList(count)

keep <- rowSums(y$counts)>=25
y <- y[keep,keep.lib.sizes=FALSE]

png("output_EdgeR_analysis/MDS.png")
plotMDS(y)
dev.off()

y <- calcNormFactors(y)

design <- model.matrix(~mESC_passage+condition,design_table)
rownames(design) <- colnames(y)

y <- estimateDisp(y,design,robust=TRUE)
png("output_EdgeR_analysis/BCV.png")
plotBCV(y)
dev.off()

fit <- glmQLFit(y, design)
png("output_EdgeR_analysis/QLDisp.png")
plotQLDisp(fit)
dev.off()

qlf <- glmQLFTest(fit,coef=4)
png("output_EdgeR_analysis/MA.png")
plotMD(qlf)
dev.off()

res = as.data.frame(topTags(qlf,n = dim(y$counts)[1]))
write.table(res,"output_EdgeR_analysis/results_EdgeR.tsv", sep="\t", quote = F, row.names = T)
FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR < 1)
