### DEXSeq ###

library(BiocParallel)

setwd("/data/manke/group/ferrari/PhD_project/reference_datasets/Franz_Vogel_dataset/RNA-Seq/output_dorsalTelencephalon_RNA-Seq_snakepipes/DEXSeq")
source("load_SubreadOutput.R")


### create info-table ###
sampleTable = data.frame(
  row.names = c("Vogel_mDorsalTelencefalon_CtrlE145_RNASeq_Rep1",
                         "Vogel_mDorsalTelencefalon_CtrlE145_RNASeq_Rep2",
                         "Vogel_mDorsalTelencefalon_CtrlE145_RNASeq_Rep3",
                         "Vogel_mDorsalTelencefalon_CtrlE145_RNASeq_Rep4",
                         "Vogel_mDorsalTelencefalon_CtrlE145_RNASeq_Rep5",
                         "Vogel_mDorsalTelencefalon_Dot1lcKOE145_RNASeq_Rep1",
                         "Vogel_mDorsalTelencefalon_Dot1lcKOE145_RNASeq_Rep2",
                         "Vogel_mDorsalTelencefalon_Dot1lcKOE145_RNASeq_Rep3",
                         "Vogel_mDorsalTelencefalon_Dot1lcKOE145_RNASeq_Rep4",
                         "Vogel_mDorsalTelencefalon_Dot1lcKOE145_RNASeq_Rep5"),
  condition = c(rep("Ctr",5), rep("Dot1lcKO",5)))

dxd <- DEXSeqDataSetFromFeatureCounts("dorsalTelencephalon_RNA-Seq.fcount.txt",
                                         flattenedfile = "gencode.vM18.annotation.DEXSeq_featurecounts.gtf",sampleData = sampleTable)





suppressPackageStartupMessages(library("DEXSeq"))

BPPARAM = MulticoreParam(workers=6)

### Filter low counts isoforms ###
dxd_filtered = dxd[apply(assay(dxd)[,1:10],1,sum) > 0,]

### Normalization ###
dxd_filtered = estimateSizeFactors(dxd_filtered)

### Dispersion Estimation ###
dxd_filtered = estimateDispersions(dxd_filtered, BPPARAM=BPPARAM)

png("output_DEXSeq/plot_Dispersion.png")
plotDispEsts(dxd_filtered)
dev.off()

### test for DEU ###
dxd_filtered = testForDEU(dxd_filtered, BPPARAM=BPPARAM)

### Estimate exon fold change ### 
dxd_filtered = estimateExonFoldChanges(dxd_filtered, fitExpToVar="condition", BPPARAM=BPPARAM)

### Summarize Results ###
dxr1 = DEXSeqResults(dxd_filtered)
# column description
mcols(dxr1)$description
# how many exons are significant with FDR < FDR_threshold ?

FDR_threshold = 0.1

table(dxr1$padj < FDR_threshold)
# how many genes are affected ?
table(tapply(dxr1$padj < FDR_threshold, dxr1$groupID, any))

# MAplot
png("output_DEXSeq/MAplot.png")
plotMA(dxr1, cex=0.8)
dev.off()

# get targets 
disreg_genes = unique(as.vector(dxr1[as.vector(dxr1$padj < FDR_threshold) %in% c(TRUE),"groupID"]))

### VISUALIZATION ### 
for ( i in disreg_genes){
  print(i)
  pdf(paste("output_DEXSeq/",unlist(strsplit(i,"[+]"))[1],"_transcripts.pdf",sep=""), width = 12, height=8)
  plotDEXSeq(dxr1, i, legend=TRUE, displayTranscripts=F, cex.axis=1.2, cex=1.3,
             lwd=2)
  dev.off()
}


### export results ###
library(rtracklayer)

anno18 <- import.gff2("/data/manke/group/ferrari/my_repository/annotations_gencode/mouse/M18/gencode.vM18.annotation.sorted.gtf", feature.type = "gene")
geneNames = data.frame(gene_ID = anno18$gene_id,
                       gene_name = anno18$gene_name)
final_gene_vector = c()
for (i in disreg_genes){
  vect_prov = unlist(strsplit(i,"[+]"))
  if (length(final_gene_vector) == 1){
    final_gene_vector = c(final_gene_vector,i)
  } else {
    for (k in vect_prov){
      final_gene_vector = c(final_gene_vector, k)
    }
  } 
}
selected = geneNames[geneNames$gene_ID %in% final_gene_vector,]

write.table(selected,"output_DEXSeq/DEXSeq_DEUgenes.tsv", quote = F, sep = "\t",row.names = F)
write.table(geneNames,"output_DEXSeq/gencode_M18_allgenes.tsv", quote = F, sep = "\t", row.names = F)
