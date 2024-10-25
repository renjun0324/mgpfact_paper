
################################################################################
#
#                                CD8_norm_count
#
################################################################################

## 1. 和GSE99254的过滤标准一样
load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE108989/1_data/norm_count.rda")
load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE108989/1_data/cellinfo.rda")
all.equal(rownames(cellinfo), colnames(norm_count))

## 2. 提取出CD8细胞
cellinfo %>% 
  filter(T_cell_type=="CD8") %>%
  dplyr::select(Cell_name) %>% 
  unlist -> cd8
cd8_count = norm_count[,cd8]
cd8_cellinfo = cellinfo[cd8,]
all.equal(rownames(cd8_cellinfo), colnames(cd8_count))

## 3. 基因名称改为symbol
load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE108989/1_data/geneinfo.rda")
rownames(geneinfo) = geneinfo$geneID
x = geneinfo[rownames(cd8_count),]
all.equal(rownames(x),rownames(cd8_count))
x = x[which(!duplicated(x$symbol)),]
x = x[!is.na(x$symbol),]

cd8_count2 = cd8_count[rownames(x),]
rownames(cd8_count2) = x$symbol

## 4. seurat流程
seurat <- CreateSeuratObject(cd8_count2, min.cells = 0,  min.features = 0, meta.data = cd8_cellinfo)
seurat <- NormalizeData(seurat, scale.factor = 100000)
seurat <- FindVariableFeatures(object = seurat, 
                          selection.method = "mean.var.plot",  
                          mean.cutoff = c(0.0125, 8),
                          dispersion.cutoff=c(0.5, Inf),
                          verbos = FALSE)
seurat <- ScaleData(object = seurat, 
               features = VariableFeatures(seurat), 
               vars.to.regress = c("Patient"), 
               do.scale = TRUE,
               do.center = TRUE,
               display.progress = TRUE)

## 5. 取出scale data去除MAIT
cd8_cellinfo %>% 
  filter(!grepl("SLC4A10",Cluster)) %>%
  dplyr::select(Cell_name) %>% 
  unlist -> s

## 6. 准备数据
cellinfo <- cd8_cellinfo[s,]
x = cellinfo$Cell_name
i = grep("P07", x)
y = x[i]
y = gsub("_", "-", y)
y2 = sapply(y, function(x) paste0(strsplit(x,"")[[1]][7:9], collapse="") )
cellinfo$position[i] = y2

## 7. save
save(cellinfo, data_matrix, seurat, file = "dat_cellinfo.rda")
save(seurat, file = "seurat.rda")

################################################################################
#
#                                CD8_norm_norm
#
################################################################################

load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE108989/1_data/norm_norm.rda")

## 基因名称改为symbol
load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE108989/1_data/geneinfo.rda")
rownames(geneinfo) = geneinfo$geneID
x = geneinfo[rownames(norm_norm),]
all.equal(rownames(x),rownames(norm_norm))
x = x[which(!duplicated(x$symbol)),]
x = x[!is.na(x$symbol),]

norm_norm2 = norm_norm[rownames(x),]
rownames(norm_norm2) = x$symbol

## 保留cd8, meta是根据既有的celltrek对象保留下来的
norm_cd8 = norm_norm2[,rownames(meta)]
data_matrix = t(norm_cd8)
cellinfo = meta[,1:25]

################################################################################
#
#                               构建CT
#
################################################################################

ct <- CreateCellTrekObject(data_matrix = data_matrix, 
                           MetaData = cellinfo, 
                           datasetTitle = "GSE108989_CD8_rmMAIT_scaled",
                           dir = "/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE108989/3_result/cd8_rmMAIT_norm_")

ct <- RUNPCA(ct, center = FALSE, scale = FALSE) ## 用scaled-data，center不需要TRUE
ct <- RUNDM(ct)
ct <- RUNtSNE(ct, npc = 1:30)
ct <- RUNUMAP(ct, npc = 1:30)
save(ct, file = "ct.rda")

## add new label
x = ct$Cluster
x[which(x=="CD8_C01-LEF1")] = "C01 T-N"
x[which(x=="CD8_C02-GPR183")] = "C02 T-CM"
x[which(x=="CD8_C03-CX3CR1")] = "C03 T-EMRA"
x[which(x=="CD8_C04-GZMK")] = "C04 T-EM"
x[which(x=="CD8_C05-CD6")] = "C05 T-RM"
x[which(x=="CD8_C06-CD160")] = "C06 IELs"
x[which(x=="CD8_C07-LAYN")] = "C07 T-EX"
x[which(x=="CD8_C08-SLC4A10")] = "C08 MAIT"
ct$celltype = x

x = ct$Cluster
x[which(x=="CD8_C01-LEF1")] = "C01 Naive CD8+ T"
x[which(x=="CD8_C02-GPR183")] = "C02 Central memory CD8+ T"
x[which(x=="CD8_C03-CX3CR1")] = "C03 EMRA CD8+ T"
x[which(x=="CD8_C04-GZMK")] = "C04 Effector memory CD8+ T"
x[which(x=="CD8_C05-CD6")] = "C05 Tissue-resident memory CD8+ T"
x[which(x=="CD8_C06-CD160")] = "C06 Intraepithelial Lymphocytes"
x[which(x=="CD8_C07-LAYN")] = "C07 Exhausted CD8+ T"
x[which(x=="CD8_C08-SLC4A10")] = "C08 Mucosal-associated invariant T"
ct$celltype2 = x

x = ct$Cluster
x[which(x=="CD8_C01-LEF1")] = "Naive CD8+ T"
x[which(x=="CD8_C02-GPR183")] = "Central memory CD8+ T"
x[which(x=="CD8_C03-CX3CR1")] = "EMRA CD8+ T"
x[which(x=="CD8_C04-GZMK")] = "Effector memory CD8+ T"
x[which(x=="CD8_C05-CD6")] = "Tissue-resident memory CD8+ T"
x[which(x=="CD8_C06-CD160")] = "Intraepithelial Lymphocytes"
x[which(x=="CD8_C07-LAYN")] = "Exhausted CD8+ T"
x[which(x=="CD8_C08-SLC4A10")] = "Mucosal-associated invariant T"
ct$celltype3 = x

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

# x = ct$celltype3
# x = factor(x, levels = d2$gse108989$celltype3)
# levels(x) = variable_lab[levels(x)]
# ct$celltype_latex = x