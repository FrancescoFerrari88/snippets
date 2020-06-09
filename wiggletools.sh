#compute AUC on regions of the bw file that are greater than 0
wiggletools apply AUC ~/ferrari/my_repository/annotations_gencode/mouse/M17/DownstreamTSS2000_ProteinCoding_limit.sorted.bed overlaps gt 0 Chronis_mESC_H3K79me2_ChIP-Seq.filtered.subtract.Input.bw Chronis_mESC_H3K79me2_ChIP-Seq.filtered.subtract.Input.bw 

#with apply_paste
wiggletools apply_paste H3K79me2_AUC_mESC_proteinCoding.txt AUC ~/ferrari/my_repository/annotations_gencode/mouse/M17/DownstreamTSS2000_ProteinCoding_limit.sorted.bed overlaps gt 0 /data/manke/group/ferrari/PhD_project/reference_datasets/mESC_Epigenome_Chronis_dataset/output_DNA-mapping_snakepipe/deepTools_ChIP/bamCompare/Chronis_mESC_H3K79me2_ChIP-Seq.filtered.subtract.Input.bw /data/manke/group/ferrari/PhD_project/reference_datasets/mESC_Epigenome_Chronis_dataset/output_DNA-mapping_snakepipe/deepTools_ChIP/bamCompare/Chronis_mESC_H3K79me2_ChIP-Seq.filtered.subtract.Input.bw
