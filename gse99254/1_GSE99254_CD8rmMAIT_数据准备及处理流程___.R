
################################################################################
#
#                                CD8_norm
#
################################################################################

## 1. 只保留CD8_norm的中心化之后的数据进行1500基因的筛选
load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE99254/3_prepare/hvg_all/cd8/norm_cd8.rda")
load("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE99254/3_prepare/hvg_all/cd8/cellinfo_cd8.rda")
all.equal(rownames(cellinfo_cd8), colnames(norm_cd8))

norm = norm_cd8
genesd <- matrixStats::rowSds(as.matrix(norm))
names(genesd) <- rownames(norm)
hvg <- names(genesd)[order(genesd,decreasing = TRUE)][1:1500]
data_matrix = norm[hvg,] %>% t

## 2. 去除MAIT
cellinfo_cd8 %>% 
  filter(!grepl("SLC4A10",celltype)) %>%
  dplyr::select(cellnames) %>% 
  unlist -> s

## 3. 准备数据
cellinfo <- cellinfo_cd8[s,]
data_matrix <- norm[hvg,s] %>% t
all.equal(rownames(data_matrix), rownames(cellinfo))

################################################################################
#
#                                CD8_norm
#
################################################################################

## 保留cd8, meta是根据既有的celltrek对象保留下来的
norm_cd8 = norm_cd8[,rownames(meta)]
data_matrix = t(norm_cd8)
cellinfo = meta[,c(1:3,11:19)]

################################################################################
#
#                                   构建CT
#
################################################################################

ct <- CreateCellTrekObject(data_matrix = data_matrix, 
                           MetaData = cellinfo, 
                           datasetTitle = "GSE99254_CD8_rmMAIT_norm",
                           dir = "/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE99254/4_result/hvg_sd1500/cd8_rmMAIT_allgene_")

ct <- RUNPCA(ct, center = FALSE, scale = FALSE) ## 用scaled-data，center不需要TRUE
ct <- RUNDM(ct)
ct <- RUNtSNE(ct, npc = 1:30)
ct <- RUNUMAP(ct, npc = 1:30)
save(ct, file = "ct.rda")

################################################################################
#
#                                  比例图
#
################################################################################

library(ggalluvial)
library(RColorBrewer)

meta = ct@MetaData
group = "Type" 
d2_label = "majortype"
getPalette = colorRampPalette(pal_d3("category20")(20))(20)
do.call(rbind, 
        lapply(c("NTC", "PTC", "TTC"), function(x){
          tmp = meta[which(meta[,group]==x),]
          if(nrow(tmp)==0){
            cell_prop = NULL
          }else{
            cell_prop <- prop.table(table(tmp[,d2_label], tmp[,group]))
            cell_prop <- data.frame(cluster = rownames(cell_prop),
                                    group = x,
                                    proportion = cell_prop[,1])
          }
          return(cell_prop)
        }) ) %>% data.frame -> cell_prop
cell_prop$group = factor(cell_prop$group, levels = c("PTC","NTC","TTC"))
cell_prop$group = as.numeric(cell_prop$group)
ggplot(cell_prop,aes(x = group, y = proportion, fill = cluster, alluvium = cluster))+
  geom_bar(stat="identity",width=0.5) +
  geom_alluvium()+
  labs(x = "", y = "Proporation", fill = "Cluster") +
  scale_fill_uchicago("default", limits = d2_levels, breaks = d2_levels, labels = d2_levels) +
  guides(fill=guide_legend(override.aes = list(size=1),ncol=1)) + 
  # guides(fill="none") + 
  rj.ftheme.small -> p
ggsave(paste0("type_majortype_proportion.pdf"), p, 
       width = 4.5 + max(strwidth(d2_levels, "inches")), height = 3.5, units = "in")

################################################################################
#
#                               增加细胞类型
#
################################################################################

x = ct$celltype
x[which(x=="CD8_C1-LEF1")] = "CD8-LEF1"
x[which(x=="CD8_C2-CD28")] = "CD8-CD28"
x[which(x=="CD8_C3-CX3CR1")] = "CD8-CX3CR1"
x[which(x=="CD8_C4-GZMK")] = "CD8-GZMK"
x[which(x=="CD8_C5-ZNF683")] = "CD8-ZNF683"
x[which(x=="CD8_C6-LAYN")] = "CD8-LAYN"
ct@MetaData$celltype2 = x