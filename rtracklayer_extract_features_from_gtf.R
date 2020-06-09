library(rtracklayer)

anno17 <- import.gff("/data/manke/group/ferrari/my_repository/annotations_gencode/mouse/M17/gencode.vM17.annotation_def.gtf", feature.type = "gene")
prot_cod_M17gencode = subset(anno17, gene_type == "protein_coding")
prot_cod_df_17 = as.data.frame(prot_cod_M17gencode)
prot_cod_df_bed_17 = prot_cod_df_17[,c(1,2,3,10,13,5,12)]

write.table(prot_cod_df_bed_17, "/data/manke/group/ferrari/my_repository/annotations_gencode/mouse/M17/protein_coding_genes_M17.bed",sep="\t", quote = F, row.names = F)
