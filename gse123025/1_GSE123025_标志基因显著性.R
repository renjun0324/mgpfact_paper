source("/share/data6/tmp/renjun/CellTrekResult/CellTrek/merge_MURP_PCA/prepare.R")
cor_theme = theme(panel.background = element_rect(fill='transparent', color="transparent"),
                  panel.border = element_rect(fill="transparent", color="transparent"),
                  panel.grid.minor=element_blank(),
                  panel.grid.major=element_blank(),
                  plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                  plot.title = element_text(size = 14, hjust = 0.5, face = "bold", colour = "black"),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  legend.key = element_rect(fill = "white"),
                  legend.title = element_text(size = 14, colour = "black"),
                  legend.text = element_text(size = 14, colour = "black"))

library(ggridges)
library(ggplotify)
library(ggnewscale)

#------------------------------------------------------------------------------
#
#                             parameters setting
#
#------------------------------------------------------------------------------

## 0. load data
load("2_pseudotime/2.4_cov_matrix/plot_gene_df.rda")
load("2_pseudotime/2.4_cov_matrix/ssdf.rda")
load("2_pseudotime/2.4_cov_matrix/curve_gene_expr_list.rda")
load("2_pseudotime/sdf.rda")
load("2_pseudotime/orig_sdf.rda")
load(paste0("5_combine_plot/",list.files("5_combine_plot","df.rda")))
load(paste0("5_combine_plot/",list.files("5_combine_plot","gall.rda")))

## 0. load tf
# scenic_path = paste0("4_differential_genes/SCENIC_500bp_lm_fdr_ct_0.05");ct_method = "lm"
scenic_path = paste0("4_differential_genes/SCENIC_500bp_aov_c_fdr_c_0.05_murpAUC_TRUE");ct_method = "aov_c"
# scenic_path = paste0("4_differential_genes/SCENIC_500bp_aov_c3&2_fdr_c_0.05");ct_method = "aov_c2"
load("2_pseudotime/2.4_cov_matrix/curve_gene_expr_list.rda")
regulon = GetRegulon(ct, scenic_path = scenic_path)
# regulon_more = GetRegulonMore(ct, scenic_path = scenic_path)

## 0. settings
# ct_method = "lm"; sig_method="fdr"; sig_label = "ct"; sig_th = 0.05
ct_method = "aov_c"; sig_method="fdr"; sig_label = "c"; sig_th = 0.05
# ct_method = "aov_c2"; sig_method="fdr"; sig_label = "c"; sig_th = 0.05
base_name = paste0(ct_method,"_",sig_method,"_",sig_label,"_",sig_th)

## 0. celltype & levels
# d2_label = "majortype"
# d2_levels = c("Naive CD8+ T", "Cytotoxic CD8+ T", "Pre Exhausted CD8+ T", 
# "Exhausted CD8+ T", "Mucosal-associated invariant T")
d2_label = "p7_ident"
d2_levels = c("P7_C0","P7_C1","P60","P7-Gpnmb+Clec7a+")
# d2_label = "celltype"
# d2_levels = c("Naive CD8+ T", "CD8 low T", "Cytotoxic CD8+ T", "Exhausted CD8+ T")
# d2_levels = c("Naive CD8+ T", "Cytotoxic CD8+ T", "Pre Exhausted CD8+ T", "Exhausted CD8+ T")
# d2_levels = c("CD8_C1-LEF1", "CD8_C2-CD28", "CD8_C3-CX3CR1", "CD8_C4-GZMK", 
# "CD8_C5-ZNF683","CD8_C6-LAYN", "CD8_C7-SLC4A10")

#-------------------------------------------------------------------------------
#                                                                              
#                              check benchmark
#                                                                              
#-------------------------------------------------------------------------------

# 0. 文章重点提到的P7-C0和P7-C1之间的差异基因
genes = c("P2ry12", "Selplg", "Tmem119", "P2ry13", "Ccr5", "Tgfbr1", "Siglech",
          "Sall1", "Spp1", "Gpnmb", "Rps14", "Igf1", "Lyz2", "Cd9", "Csf1", "Clec7a","Lpl",
          "Cd63", "Ctsd", "Apoe", "Lgals3", "Ctsl", "Lilrb4", "Fth1", "Fabp5", "Ank",
          "Timp2", "Ctsz", "Itgax")
# genes = c("Clec7a", "Gpnmb", "Igf1", "Tmem119", "P2ry12", "Tgfbr1")

# 1. 导入文章结果显著的基因
lapply(c(2,3,5), function(i){
  x = xlsx::read.xlsx("../../0_paper/mmc6.xlsx", sheetIndex = i)
  x$gene[which(x$p_val_adj<0.05)]
}) -> genelist
unique(unlist(genelist)) -> genes

# 2. 和现有的结果取交集查看效果
for(i in 1:L){
  y = GetSigGene(ct@BeautGene, i, ct_method, sig_method, sig_label, sig_th)
  x = intersect(y, genes)
  cat(x, "\n")
  cat(length(x), "\n")
}

lapply(1:3, function(i){
  GetSigGene(ct@BeautGene, i, ct_method, sig_method, sig_label, sig_th)
  }) -> genelist