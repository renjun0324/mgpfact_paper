

#------------------------------------------------------------------------------
#
#                             original paper
#
#------------------------------------------------------------------------------

# The Seurat package was used to perform unsupervised clustering analysis on scRNA-seq data (Macosko et al., 2015). 
# Briefly, gene counts for cells that passed QC were normalized to the total expression and log-transformed, and then highly variable genes were
# detected (y.cutoff = 0.5). 
# Depending on the analysis, cell cycle effect could be regressed out using the ScaleData function. Using
# highly variable genes as input, principal component analysis was performed on the scaled data in order to reduce dimensionality.
# Statistically significant principal components were determined by using the JackStrawPlot function. These principal components
# were used to compute the distance metric, which then generated cell clusters. Non-linear dimensional reduction (tSNE) was used
# to visualize clustering results. Differentially expressed genes were found using the FindAllMarkers (or FindMarkers) function that
# ran Wilcoxon rank sum tests.

#------------------------------------------------------------------------------
#
#                                数据读取
#
#------------------------------------------------------------------------------

require(data.table)
require(xlsx)
require(dplyr)


count <- fread("1_data/GSE123025_Single_myeloid_1922_cells_processed_data.csv.gz",data.table=FALSE)
rownames(count) <- count[,1]
count <- count[,-1]

cellinfo <- read.xlsx("1_data/GSE123025_meta_data_scRNA_Li.xlsx", sheetIndex = 1)
rownames(cellinfo) = cellinfo[,1]
cellinfo <- cellinfo[,-1]
cellinfo2 <- cellinfo[colnames(count),]
eql(colnames(count), rownames(cellinfo2))

cellinfo2$day = "a"
cellinfo2$day[which(cellinfo2$gate %in% c("cKit-CD45+"))] = "E14.5"
cellinfo2$day[which(cellinfo2$gate %in% c("CD45hi", "CD45low"))] = "P7"
cellinfo2$day[which(cellinfo2$gate %in% c("Tmem119-", "Tmem119+"))] = "P60"
cellinfo3 = cellinfo2[which(cellinfo2$ori_cluster_ident!="NA"),]

cellinfo = cellinfo3
count = count[,rownames(cellinfo3)]
x = cellinfo$ori_cluster_ident
x[which(x=="0")] = "0 Homeostatic microglia"
x[which(x=="1")] = "1 Postnatal PAM"
x[which(x=="2")] = "2 Postnatal immature microglia"
x[which(x=="3")] = "3 G2/M phase microglia"
x[which(x=="4")] = "4 S phase microglia"
x[which(x=="5")] = "5 Embryonic-like microglia"
x[which(x=="6")] = "6 IEG+ microglia"
x[which(x=="7")] = "7 Adult choroid plexus Mφ"
x[which(x=="8")] = "8 Patrolling monocytes"
x[which(x=="9")] = "9 Postnatal choroid plexus Mφ"
x[which(x=="10")] = "10 CD209a+ monocytes"
x[which(x=="11")] = "11 Inflammatory monocytes"
x[which(x=="12")] = "12 Neutrophils"
x[which(x=="13")] = "13 Choroid plexus epithelial cells"
x[which(x=="14")] = "14 NK cells"
cellinfo$celltype = x 

save(cellinfo, file = "1_data/cellinfo.rda")
save(count, file = "1_data/count.rda")

#------------------------------------------------------------------------------
#
#                                数据读取
#
#------------------------------------------------------------------------------

## seurat
seurat <- CreateSeuratObject(count, 
                                     project = "s", 
                                     min.cells = 0, 
                                     min.features = 0, 
                                     assay = "RNA",
                                     names.field = 1,
                                     names.delim = "_",
                                     meta.data = cellinfo)

seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(object = seurat, 
                                       selection.method = "mean.var.plot",  
                                       # mean.cutoff = c(0.0125, 3),
                                       dispersion.cutoff=c(0.5, Inf),
                                       verbos = FALSE)
seurat <- ScaleData(seurat)
save(seurat, file = "2_prepare/seurat.rda")

## save
sdata <- seurat@assays$RNA@data[VariableFeatures(seurat),] 
sdata <- t(as.matrix(sdata))
cellinfo <- data.frame(row.names = rownames(seurat@meta.data),
                       cell_name = rownames(seurat@meta.data),
                       seurat@meta.data,
                       truetime = seurat@meta.data$day,
                       stringsAsFactors = FALSE)

## rename truetime
x = do.call(rbind,strsplit(cellinfo$truetime, "P"))[,2]
x = do.call(rbind,strsplit(x, "E"))[,2]
cellinfo$truetime = as.numeric(x)

save(sdata, file = "2_prepare/sdata.rda")
save(cellinfo, file = "2_prepare/cellinfo.rda")

#------------------------------------------------------------------------------
#
#                             按照文章提取部分细胞
#
#------------------------------------------------------------------------------

load("2_prepare/seurat_1816.rda")

# 1. choose
p7_microg = fread("1_data/GSE123024_Gpnmb_Clec7a_143_Microglia_processed_data.csv.gz", data.table = FALSE)
rownames(p7_microg) = p7_microg[,1]
p7_microg = p7_microg[,-1]

load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE123025/1_data/count.rda")
load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE123025/2_prepare/cellinfo.rda")
tmp = sdata[which(cellinfo$P60_MG_cluster_ident=="homeostatic"),"Clec7a"]
p60_ind = which(tmp==0) %>% names
p7_ind = which(cellinfo$P7_MG_cluster_ident %in% c("P7_C0", "P7_C1"))

# 2. combine
count = do.call(cbind, list(count[,p60_ind], count[,p7_ind], p7_microg))
cellinfo = data.frame(row.names = colnames(count),
                      day = c(rep("P60", length(p60_ind)),
                              rep("P7", length(p7_ind)),
                              rep("P7-Gpnmb+Clec7a+", ncol(p7_microg))) )

# 3. rename truetime
x = do.call(rbind,strsplit(cellinfo$day,"-"))[,1]
cellinfo$truetime = do.call(rbind,strsplit(x,"P"))[,2]%>% as.numeric

# 4. add other cellinfo
ind = intersect(colnames(seurat), rownames(cellinfo))
cellinfo$celltype = "P7-Gpnmb+Clec7a+"
cellinfo[ind,"celltype"] = seurat@meta.data[ind,"celltype"]
cellinfo$p7_ident = "P7-Gpnmb+Clec7a+"
p7_names = rownames(seurat@meta.data)[which(seurat$P7_MG_cluster_ident %in% c("P7_C0", "P7_C1"))]
p7_names = intersect(ind, p7_names)
cellinfo[p7_names,"p7_ident"] = seurat@meta.data[p7_names,"P7_MG_cluster_ident"]
cellinfo[which(cellinfo$day=="P60"),"p7_ident"] = "P60"

# 5. save
cellinfo = cellinfo[colnames(count),]
all.equal(colnames(count), rownames(cellinfo))
save(count, cellinfo, file = "2_prepare/paper_for_monocle.rda")

# 6. seurat
seurat <- CreateSeuratObject(count, 
                             project = "s", 
                             min.cells = 0, 
                             min.features = 0, 
                             assay = "RNA",
                             names.field = 1,
                             names.delim = "_",
                             meta.data = cellinfo)

seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(object = seurat, 
                               selection.method = "mean.var.plot",  
                               # mean.cutoff = c(0.0125, 3),
                               dispersion.cutoff=c(0.5, Inf),
                               verbos = FALSE)
seurat <- ScaleData(seurat)
save(seurat, file = "2_prepare/seurat_paper_for_monocle.rda")

## 7. save
sdata <- seurat@assays$RNA@data[VariableFeatures(seurat),] 
sdata <- t(as.matrix(sdata))
cellinfo <- data.frame(row.names = rownames(seurat@meta.data),
                       cell_name = rownames(seurat@meta.data),
                       seurat@meta.data,
                       stringsAsFactors = FALSE)
save(sdata, file = "2_prepare/paper_for_monocle/sdata.rda")
save(cellinfo, file = "2_prepare/paper_for_monocle/cellinfo.rda")
