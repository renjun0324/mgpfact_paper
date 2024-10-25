
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
# load("2_pseudotime/2.4_cov_matrix/plot_gene_df.rda")
# load("2_pseudotime/2.4_cov_matrix/ssdf.rda")
# load("2_pseudotime/2.4_cov_matrix/curve_gene_expr_list.rda")
# load("2_pseudotime/sdf.rda")
# load("2_pseudotime/orig_sdf.rda")
# load(paste0("5_combine_plot/",list.files("5_combine_plot","df.rda")))
# load(paste0("5_combine_plot/",list.files("5_combine_plot","gall.rda")))
ssdf = ct@GPR$curve_sdf
curve_gene_expr_list = ct@GPR$gene_expr
sdf = GetMURPInfo(ct)

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

genes = c("P2ry12", "Selplg", "Tmem119", "P2ry12", "P2ry13", "Ccr5", "Tgfbr1", "Siglech",
          "Sall1", "Spp1", "Gpnmb", "Rps14", "Igf1", "Lyz2", "Cd9", "Csf1", "Clec7a","Lpl",
          "Cd63", "Cd63", "Ctsd", "Apoe", "Lgals3", "Ctsl", "Lilrb4", "Fth1", "Fabp5", "Ank",
          "Timp2", "Ctsz", "Itgax")

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
  cat(length(x), "\n")
}

# #-------------------------------------------------------------------------------
# #                                                                              
# #                           expression + GRP ~ C
# #                                                                              
# #-------------------------------------------------------------------------------
# 
# # genes = c("Stat1","Clec7a","Hif1a","Cd28","Lmna","Abca1","Igf1","Lyz2","Plek","Ldha")
# genes = c("Tmem119", "Selplg", "P2ry12", "Tgfbr1")
# for(l in 1:3){
#   
#   tb = sdf[1,paste0("Tb_",l)]
#   
#   lapply(genes, function(g){
#     cat(g, "\n")
#     df = data.frame(sdf, gene = ct@MURP$Recommended_K_cl$centers[,g])
#     gpr_df <- data.frame(plot_gene_df[,paste0("L",l,"_",g),drop=FALSE],
#                          T = ssdf$T,
#                          C = as.factor(ssdf[,paste0("C_",l)]) ) %>% melt(id = c("T","C"))
#     gpr_df$variable <- stringr::str_split_fixed(gpr_df$variable, "_", 2)[,2]
#     gpr_df1 <- gpr_df[which(gpr_df$C==1),]
#     gpr_df2 <- gpr_df[which(gpr_df$C==2),]
#     
#     # plot
#     df[,paste0("C0_",l)] = factor(df[,paste0("C0_",l)])
#     ggplot() + 
#       geom_vline(xintercept = tb, colour = "#990000", linetype = "dashed") +
#       geom_point(data = df, size = 3, alpha = 0.5, shape = 21,
#                  aes(x = T, y = gene, fill =  p7_ident), colour = "white" ) +
#       geom_smooth(data = gpr_df, aes(x = T, y = value, colour = C),
#                   alpha = 0.7, size = 1, method = loess, formula = y ~ x, se = FALSE, na.rm = TRUE) +
#       # scale_fill_d3("category20") +
#       # scale_color_d3("category20") +
#       scale_fill_manual(values = col_values, limits = d2_levels, breaks = d2_levels, labels = d2_levels, drop=FALSE ) +
#       scale_colour_manual(values = col_values, limits = d2_levels, breaks = d2_levels, labels = d2_levels, drop=FALSE ) +
#       
#       scale_color_manual(values = c(up = "#374e55ff", down = "#df8f44ff"),
#                          # values = c(up = "#159164", down = "#c26c35"),
#                          # values = c(up = "#D074C1", down = "#84d7e1ff"),
#                          breaks = c("up", "down"),
#                          labels = c("Up in\nPhase 0->1\n", "Up in\nPhase 0->2\n"),
#                          guide = guide_legend(override.aes = list(width = 4))) +
#       
#       labs(title = g, x = "Pseudotime", y = "Expression", color = "Bif") +
#       guides(shape = "none", colour = "none") + 
#       rj.ftheme
#     # ggplot() + 
#       # geom_point(data = df1, aes(x = UMAP1, y = UMAP2), 
#       #            colour = "grey80", alpha = 0.2, size = 2.3) + 
#       # geom_point(data = df2, aes(x = UMAP1, y = UMAP2, colour = celltype), 
#       #            alpha = 0.6, size = 2.3) + 
#       # scale_colour_d3("category20", breaks = d2_levels, labels = d2_levels,drop=FALSE ) +
#      
#       guides(colour = guide_legend(override.aes=list(size=3,alpha=0.8))) +
#       labs(colour = NULL, x = "UMAP-1", y = "UMAP-2") +
#       rj.ftheme -> p_murp_umap_celltype
#     
#   }) -> plist
#   p = patchwork::wrap_plots(plist, ncol = 4)
#   
#   ## save
#   if(length(genes)==1){
#     ggsave(paste0("4_differential_genes/1_exprT_GPRsmooth_L", l, "_", paste(genes, collapse = "_"),".pdf"), p, 
#            width = 11, height = length(genes)*4, dpi = 150)
#   }else{
#     ggsave(paste0("4_differential_genes/1_exprT_GPRsmooth_L",l,"_", paste(genes, collapse = "_"),".pdf"), p, 
#            width = 11, height = length(genes)*1, dpi = 150)
#   }
# }
# 
# #-------------------------------------------------------------------------------
# #                                                                              
# #                        expression + GRP ~ Celltype
# #                                                                              
# #-------------------------------------------------------------------------------
# 
# # genes = c("Stat1","Clec7a","Hif1a","Cd28","Lmna","Abca1","Igf1","Lyz2","Plek","Ldha")
# # genes = c("Tmem119, Selplg, P2ry12, Tgfbr1")
# for(l in 1:3){
#   
#   tb = sdf[1,paste0("Tb_",l)]
#   
#   ## 转录因子 
#   # genes = regulon$tf[[l]]
#   
#   ## density
#   lapply(genes, function(g){
#     cat(g, "\n")
#     
#     # expr
#     # df = data.frame(orig_sdf, gene = ct@assay$data_matrix[,g])
#     df = data.frame(sdf, gene = ct@MURP$Recommended_K_cl$centers[,g])
#     
#     # GRP
#     gpr_df <- data.frame(plot_gene_df[,paste0("L",l,"_",g),drop=FALSE],
#                          T = ssdf$T,
#                          C = as.factor(ssdf[,paste0("C_",l)]) ) %>% melt(id = c("T","C"))
#     gpr_df$variable <- stringr::str_split_fixed(gpr_df$variable, "_", 2)[,2]
#     gpr_df1 <- gpr_df[which(gpr_df$C==1),]
#     gpr_df2 <- gpr_df[which(gpr_df$C==2),]
#     
#     # plot
#     df[,paste0("C0_",l)] = factor(df[,paste0("C0_",l)])
#     ggplot() + 
#       geom_vline(xintercept = tb, colour = "#990000", linetype = "dashed") +
#       geom_point(data = df, 
#                  aes(x = T, y = gene, color = get(d2_label), shape = get(d2_label) ),
#                  size = 2.3, alpha = 0.7) + 
#       geom_smooth(data = gpr_df1, aes(x = T, y = value), color = "grey20",
#                   alpha = 0.5, size = 0.8, linetype = "dashed",
#                   method = loess, formula = y ~ x, se = FALSE, na.rm = TRUE) +
#       geom_smooth(data = gpr_df2, aes(x = T, y = value), color = "grey50",
#                   alpha = 0.2, size = 0.8, linetype = "dashed",
#                   method = loess, formula = y ~ x, se = FALSE, na.rm = TRUE) +
#       # geom_line(data = gpr_df1, aes(x = T, y = value, colour = C),
#       #             alpha = 0.7, size = 1) +
#       # geom_line(data = gpr_df2, aes(x = T, y = value, colour = C),
#       #             alpha = 0.7, size = 1) +
#       scale_color_d3("category20") + 
#       labs(title = g, x = "Pseudotime", y = "Expression", color = "Status", shape = "Status") +
#       guides(colour = guide_legend(override.aes=list(size=4, alpha = 0.6))) +
#       rj.ftheme
#     
#   }) -> plist
#   p = patchwork::wrap_plots(plist, ncol = 4, guides = "collect")
#   
#   ## save
#   if(length(genes)==1){
#     ggsave(paste0("4_differential_genes/1_exprT_GPRsmooth_",d2_label,"L", l, "_", paste(genes, collapse = "_"),".pdf"), p, 
#            width = 14, height = length(genes)*4, dpi = 150)
#   }else{
#     ggsave(paste0("4_differential_genes/1_exprT_GPRsmooth_",d2_label,"L", l, "_", paste(genes, collapse = "_"),".pdf"), p, 
#            width = 14, height = length(genes)*1, dpi = 150)
#   }
# }