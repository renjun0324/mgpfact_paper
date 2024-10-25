
## 确认文章用来做轨迹的基因
lapply(d2$gse108989$Cluster, function(i){
  cat(i, "\n")
  genes = xlsx::read.xlsx("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE108989/0_paper/Supplementary Table 5.xlsx", 
                          sheetName = i)
  colnames(genes) = genes[1,]
  genes = genes[-1,1:3]
  genes = genes[!is.na(genes[,2]),]
  genes[,3] = as.numeric(genes[,3])
  genes = genes[order(genes[,3], decreasing = T),]
  return(genes)
}) -> glist
gdf = do.call(rbind, glist)
gdf = gdf[order(gdf[,3], decreasing = T),]
g = unique(gdf[1:700,2])

data_matrix = ct@assay$data_matrix[,g]
metadata = ct@MetaData[,c(1:16,27:34)]

ct <- CreateCellTrekObject(data_matrix = data_matrix, 
                           MetaData = metadata,  
                           datasetTitle = "paper_monocle_degs")

x = ct$Cluster
x[which(x=="CD8_C01-LEF1")] = "CD8-LEF1"
x[which(x=="CD8_C02-GPR183")] = "CD8-GPR183"
x[which(x=="CD8_C03-CX3CR1")] = "CD8-CX3CR1"
x[which(x=="CD8_C04-GZMK")] = "CD8-GZMK"
x[which(x=="CD8_C05-CD6")] = "CD8-CD6"
x[which(x=="CD8_C06-CD160")] = "CD8-CD160"
x[which(x=="CD8_C07-LAYN")] = "CD8-LAYN"
x[which(x=="CD8_C08-SLC4A10")] = "CD8-SLC4A10"
ct$celltype4 = x

save(ct, file = "ct.rda")

data_matrix = ct@assay$data_matrix
meta_data = ct@MetaData
i = which(meta_data$celltype3!="Mucosal-associated invariant T")
data_matrix = data_matrix[i,]
meta_data = meta_data[i,]
ct <- CreateCellTrekObject(data_matrix = data_matrix, 
                           MetaData = meta_data,  
                           datasetTitle = "paper_monocle_degs")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                              添加一个新标签
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = ct$celltype4
x[which(x=="CD8-LEF1")] = "CD8-LEF1 (Naive CD8+ T)"
x[which(x=="CD8-GPR183")] = "CD8-GPR183 (Central memory CD8+ T)"
x[which(x=="CD8-CX3CR1")] = "CD8-CX3CR1 (EMRA CD8+ T)"
x[which(x=="CD8-GZMK")] = "CD8-GZMK (Effector memory memory CD8+ T)"
x[which(x=="CD8-CD6")] = "CD8-CD6 (Tissue-resident memory CD8+ T)"
x[which(x=="CD8-LAYN")] = "CD8-LAYN (Exhausted CD8+ T)"
x[which(x=="CD8-CD160")] = "CD8-CD160 (Intraepithelial Lymphocytes)"
x = factor(x, levels = d2[["gse108989"]][["celltype5"]])
ct$celltype5 = x

levels(x) = variable_lab[levels(x)]
ct@MetaData$celltype6 = x
ct = GetMURPMapLabel(ct, labels = c("celltype6"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                             增加一个tpm的仓位
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE108989/1_data/norm_tpm_cd8.rda")
load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE108989/1_data/geneinfo.rda")
rownames(geneinfo) = geneinfo$geneID
x = geneinfo[rownames(norm_tpm_cd8),]
all.equal(rownames(x),rownames(norm_tpm_cd8))
x = x[which(!duplicated(x$symbol)),]
x = x[!is.na(x$symbol),]

norm_tpm_cd8_2 = norm_tpm_cd8[rownames(x),]
rownames(norm_tpm_cd8_2) = x$symbol
x = norm_tpm_cd8_2 %>% t
ct@assay$tpm_all_gene = x[rownames(ct@assay$data_matrix_all_gene), colnames(ct@assay$data_matrix_all_gene)]
