# Biomart script

library(biomaRt)

listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)


filters = listFilters(ensembl)
filters[1:5,]

attributes = listAttributes(ensembl)
attributes[1:5,]


genes = getBM(attributes=c("chromosome_name",
                   "start_position",
                   "end_position",
                   "ensembl_gene_id_version",
                   "external_gene_name",
                   "strand",
                   "gene_biotype"), 
      filters = c("biotype","chromosome_name"), 
      values = list(c("protein_coding"),c("1")), 
      mart = ensembl)
