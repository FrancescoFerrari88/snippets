# ####################################
# ## FUNCTIONAL ENRICHMENT ANALYSIS ##
# ####################################
# 

### Foreground VS Background ###

library(clusterProfiler)
library(org.Mm.eg.db)


### load table 
res = read.table("....")

res_no.na = na.omit(res)
res_lfc_cutoff_no.na = na.omit(res_lfc_cutoff)

target_genes_noLFCthr = res_no.na$symbol[res_no.na$padj < 0.05]
target_genes_withLFCthr = res_lfc_cutoff_no.na$symbol[res_lfc_cutoff_no.na$padj < 0.05]
background_genes = res_no.na$symbol

keytypes(org.Mm.eg.db)

### without LFC cutoff ###
eg = bitr(target_genes_noLFCthr, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
universe = bitr(background_genes, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
lfc = as.vector(res_no.na$log2FoldChange[match(eg$SYMBOL, res_no.na$symbol)])
names(lfc) = eg$ENTREZID

ego <- enrichGO(gene          = eg[["ENTREZID"]],
                universe      = universe[["ENTREZID"]],
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)

pdf("output_DESeq2_analysis/enrichment_analysis_allDEgenes/dotplot_DEgenes_noLFCthr.pdf")
dotplot(ego, showCategory = 20)
dev.off()
pdf("output_DESeq2_analysis/enrichment_analysis_allDEgenes/barplot_DEgenes_noLFCthr.pdf")
barplot(ego, showCategory = 20)
dev.off()
pdf("output_DESeq2_analysis/enrichment_analysis_allDEgenes/emapplot_DEgenes_noLFCthr.pdf")
emapplot(ego, showCategory = 50)
dev.off()
pdf("output_DESeq2_analysis/enrichment_analysis_allDEgenes/cnetplot_DEgenes_noLFCthr.pdf", width = 10, height = 8)
cnetplot(ego, showCategory = 8, categorySize="pvalue", foldChange = lfc)
dev.off()
pdf("output_DESeq2_analysis/enrichment_analysis_allDEgenes/cnetplot_DEgenes_noLFCthr_20.pdf", width = 20, height = 20)
cnetplot(ego, showCategory = 20, categorySize="pvalue", foldChange = lfc)
dev.off()

### get table
as.data.frame(ego)





### Compare Clusters ###

cl = read.csv("genes_N1-N2.txt",sep="\t")
if(grepl("\\.[0-9]{1,2}",cl$name)){cl$name<-gsub("\\.[0-9]+","",cl$name)}
cl_1 = as.vector(subset(cl, deepTools_group=="cluster_1")[,"name"])
cl_2 = as.vector(subset(cl, deepTools_group=="cluster_2")[,"name"])
cl_3 = as.vector(subset(cl, deepTools_group=="cluster_3")[,"name"])
cl_4 = as.vector(subset(cl, deepTools_group=="cluster_4")[,"name"])
cl_5 = as.vector(subset(cl, deepTools_group=="cluster_5")[,"name"])

clust_tab = list(cluster_1 = cl_1,
                 cluster_2 = cl_2,
                 cluster_3 = cl_3,
                 cluster_4 = cl_4,
                 cluster_5 = cl_5)

clust_tab_def = list()
# clust_tab = c("cluster_1","cluster_2","cluster_3","cluster_4","cluster_5")
print(names(clust_tab))

for (i in 1:5){
  v = bitr(clust_tab[[i]], fromType="ENSEMBL", toType=c("ENTREZID","ENSEMBL","SYMBOL"), OrgDb="org.Mm.eg.db")
  clust_tab_def[[i]] = v$ENTREZID
  print(names(clust_tab)[i])
  print(head(clust_tab_def[[i]]))
}
names(clust_tab_def) = names(clust_tab)


ck_Tab_def <- compareCluster(geneCluster = clust_tab_def,
                             fun = "enrichGO",
                             OrgDb="org.Mm.eg.db",
                             ont="BP",
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.1,
                             pvalueCutoff = 0.1,
                             readable      = TRUE)

pdf("GO_enrichment_clusters.pdf", width = 11, height = 6)
dotplot(ck_Tab_def)
dev.off()

