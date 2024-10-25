
genes = xlsx::read.xlsx("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE99254/0_paper/41591_2018_45_MOESM5_ESM_clustergenes.xlsx", sheetIndex = 1)
idd = grep("CD8(?!.*SLC4A10)", genes[,3], perl = TRUE, value = FALSE)
g = genes[idd, 2] %>% unique

load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE99254/4_result/hvg_sd1500/cd8_rmMAIT_allgene_witho/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE99254/4_result/hvg_sd1500/cd8_rmMAIT_allgene_withoutroot/ct.rda")
g = intersect(g,colnames(ct@assay$data_matrix)) 

data_matrix = ct@assay$data_matrix[,g]
metadata = ct@MetaData[,c(1,3:12)]

ct <- CreateCellTrekObject(data_matrix = data_matrix, 
                           MetaData = metadata,  
                           datasetTitle = "GSE99254_CD8_norm",
                           dir = "paper_monocle_degs")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                              添加一个新标签
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = ct$celltype2
x[which(x=="CD8-LEF1")] = "CD8-LEF1 (Naive CD8+ T)"
x[which(x=="CD8-CD28")] = "CD8-CD28 (Naive CD8+ T)"
x[which(x=="CD8-CX3CR1")] = "CD8-CX3CR1 (Cytotoxic CD8+ T)"
x[which(x=="CD8-GZMK")] = "CD8-GZMK (Pre Exhausted CD8+ T)"
x[which(x=="CD8-ZNF683")] = "CD8-ZNF683 (Pre Exhausted CD8+ T)"
x[which(x=="CD8-LAYN")] = "CD8-LAYN (Exhausted CD8+ T)"
x = factor(x, levels = d2[["gse99254"]][["celltype5"]])
ct$celltype5 = x

levels(x) = variable_lab[levels(x)]
ct@MetaData$celltype6 = x
ct = GetMURPMapLabel(ct, labels = c("celltype5", "celltype6"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                             增加一个tpm的仓位
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE99254/3_prepare/hvg_all/cd8/tpm_cd8.rda")
x = t(tpm_cd8)
ct@assay$tpm_all_gene = x[rownames(ct@assay$data_matrix_all_gene), colnames(ct@assay$data_matrix_all_gene)]

