
## 利用weight基因做ICI的ROC曲线

# julia_home = "/public/home/renjun/tool/julia-1.6.6/bin"
# celltrek_root <- "/share/data6/tmp/renjun/CellTrekResult/CellTrek/celltrekr/R"
# invisible(sapply(list.files(celltrek_root, full.names = T), source))
# celltrek_root <- "/share/data6/tmp/renjun/CellTrekResult/CellTrek/celltrek_code_all"
# source(paste0(celltrek_root,"/prepare.R"))

library(rstatix)
library(patchwork)
library(pROC)
library(MGPfactR)
library(magrittr)
library(ggsci)
library(ggplotify)
celltrek_root <- "/share/data6/tmp/renjun/CellTrekResult/CellTrek/celltrek_code_all"
dir.create("5_combine_plot/ici")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                  准备数据
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# sdf = GetMURPInfo(ct)
# data_matrix = ct@MURP$Recommended_K_cl$centers
# sdf = ct@MetaData
# data_matrix = ct@assay$data_matrix

gene_weight = GetGeneWeight(ct)
zscore = T
load(paste0(celltrek_root, "/ici/deseq2_normalize_data.rda"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                           高权重基因的加权轨迹得分
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zscore = TRUE
# as <- c("Cho2020","Auslander2018","TCGA-SKCM") # gse99254
as <- c("Jung2019","Cho2020","Auslander2018","TCGA-SKCM") # gse99254
weight_cut = 0.05

lapply(1:3, function(l){
  gw = gene_weight[,l]
  weight_gene = names(gw)[which(abs(gw) > weight_cut)]
  
  lapply(as, function(index){
    
    indd = which(clin$author==index)
    ici_group = clin[indd, "Response"]
    gs = intersect(weight_gene, rownames(normalize_data))
    wg = gw[gs] # 包含权重和基因名称的向量
    
    ici_expr = normalize_data[gs, indd, drop=F]
    if(zscore){
      ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
    }
    for(xi in rownames(ici_expr)){
      ici_expr[xi,] = ici_expr[xi,] * wg[xi]
    }
    
    ici_means = colMeans(ici_expr)
    rocc = roc(ici_group, ici_means, levels=c("Response", "Non_response"), ci=TRUE)
    return(rocc)
    
  }) -> cl
  
  names(cl) = as
  glab = sapply(seq_along(cl), function(i){
    paste0(names(cl)[i], " (AUC=",as.numeric(auc(cl[[i]])) %>% round(4),")")
  }) %>% set_names(as)
  
  pROC::ggroc(cl, size = 1, aes=c("colour"), alpha = 1, legacy.axes = TRUE) +
    geom_abline(slope=1, intercept=0, color = "#3e3e23ff", alpha = 0.8, size = 0.6, linetype = "dashed") +
    labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", title = paste0("Trajectory ", l)) +
    scale_color_aaas(name = "", breaks = names(glab), labels = glab) +
    guides(linetype="none") +
    theme(panel.background = element_rect(fill='transparent', color="black", size = 0.3),
          strip.text = element_text(size = 10),
          strip.background = element_rect(colour = "black", fill = "transparent", size = 0.25),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.border = element_rect(fill='transparent', color='black', size = 0.2),
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
          plot.title = element_text(vjust = 0, size = 16, colour = 'black'),
          legend.key = element_rect( fill = "white"),
          legend.position = c(0.65,0.2),
          legend.background = element_rect(fill='transparent'),
          legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
          legend.key.size = unit(0.6, "cm"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(vjust = 0, size = 16, colour = 'black'),
          axis.title.y = element_text(vjust = 2, size = 18, colour = 'black',margin = margin(r = 5)),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 16, colour = 'black', margin = margin(t = 6)),
          axis.text.y = element_text(vjust = 0.5, size = 16, colour = 'black',margin = margin(r = 6)) ) -> p2
  p2 = as.ggplot(as.grob(p2))
  
  return(p2)
}) -> pl
p = wrap_plots(pl, nrow = 1) +
  plot_annotation(title = paste0("Weight CUT: ",weight_cut)) &
  theme(plot.title = element_text(size = 14, face="bold"))

pdf(paste0("5_combine_plot/ici/ici_weight_gene_roc_zscore",zscore,"_trajectory_weightmeans_0.05_3个确定的数据.pdf"), width = 15, height = 5.5)
p
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                           高权重基因的加权轨迹得分
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zscore = TRUE
as <- c("Cho2020","Auslander2018","TCGA-SKCM") # gse99254
# as <- c("Jung2019","Cho2020","Auslander2018","TCGA-SKCM") # gse99254
weight_cut = 0.05

plist = list()
for(l in 1:2){
  
  ## (1).  计算得分返回分组
  gw = gene_weight[,l]
  weight_gene = names(gw)[which(abs(gw) > weight_cut)]
  
  lapply(as, function(index){
    
    indd = which(clin$author==index)
    ici_group = clin[indd, "Response"]
    gs = intersect(weight_gene, rownames(normalize_data))
    wg = gw[gs] # 包含权重和基因名称的向量
    ici_expr = normalize_data[gs, indd, drop=F]
    if(zscore){
      ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
    }
    for(xi in rownames(ici_expr)){
      ici_expr[xi,] = ici_expr[xi,] * wg[xi]
    }
    ici_means = colMeans(ici_expr)
    df = data.frame(expr = ici_means, 
                    group = factor(ici_group,levels=c("Response", "Non_response")))
    df$author = index
    return(df)
  }) -> lx
  dff = do.call(rbind, lx) %>% as.data.frame
  
  dff$author[which(dff$author=="Jung2019")] = "NSCLC-Jung"
  dff$author[which(dff$author=="Auslander2018")] = "Melanoma-Auslander"
  dff$author[which(dff$author=="Cho2020")] = "NSCLC-Cho"
  
  ## (2).  计算显著性
  annotation_df = get_pval_df(dff = dff, 
                              attr_name = c("author"),
                              value_name = "expr", 
                              group_name = "group", 
                              test_method = "t.test")
  annotation_df$label = sapply(annotation_df$pval, function(p){
    p2 = as.numeric(p)
    if(is.na(p2)){
      p
    }else{
      pval_lab(p2, lab = T)
    }
  })
  
  ## (3).  画图
  p <- ggplot(dff) + 
    geom_boxplot(aes(x = group, y = expr, fill = group),
                 position = "dodge", width=0.4, size = 0.1,outlier.shape = NA) +
    labs(x = "", y = "Weighted mean of HWGs", fill = "") +
    scale_fill_grey() +
    facet_wrap(~author, nrow = 1) + 
    guides(fill = "none") + 
    scale_y_continuous(limits = c(-0.1,0.1)) +
    geom_signif( data = annotation_df,
                 aes(xmin = X1, xmax = X2, annotations = TeX(label, output = "character"),
                     y_position = y_position-0.03),
                 tip_length = 0,
                 size = 0.06,
                 textsize = 5,
                 vjust = -0.1,
                 linetype = "dashed",
                 test = wilcox.test,
                 manual = TRUE,
                 parse = TRUE) +
    theme(panel.background = element_rect(fill='transparent', color="black", size = 0.5),
          strip.text = element_text(size = 13),
          strip.background = element_rect(colour = "black", fill = "transparent", size = 0.5),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.border = element_rect(fill='transparent', color='black', size = 0.5),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(vjust = -1.5, size = 13, colour = 'black'),
          axis.title.y = element_text(vjust = 1.5, size = 13, colour = 'black'),
          axis.text.y.left = element_text(vjust = 0, hjust = 1, size = 13, colour = 'black'),
          axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, size = 13, colour = 'black'),
          # axis.text.x = element_text(size = 3.7,colour = 'black'),
          strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")),
          strip.text.y = element_text(margin = margin(0,0.05,0,0.05, "cm")))
  plist = append(plist, list(p))
}
p = wrap_plots(plist, nrow = 2)
ggsave("independent_traj_box.pdf", p, width = 5.5, height = 7) # gse99254
# ggsave("independent_traj_box.pdf", p, width = 7, height = 7) # gse108989

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                               计算ROC的P值
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

get_pval_roc <- function(r){
  v <- var(r)  # 计算方差
  se <- sqrt(v)  # 从方差计算标准误差
  b <- r$auc - .5
  z <- (b / se)
  p_value <- 2 * pt(-abs(z), df=Inf)
  return(p_value)
}

zscore = TRUE
as <- c("Cho2020","Auslander2018","TCGA-SKCM") # gse99254
# as <- c("Jung2019","Cho2020","Auslander2018","TCGA-SKCM") # gse108989
weight_cut = 0.05

lapply(1:3, function(l){
  gw = gene_weight[,l]
  weight_gene = names(gw)[which(abs(gw) > weight_cut)]
  
  lapply(as, function(index){
    
    indd = which(clin$author==index)
    ici_group = clin[indd, "Response"]
    gs = intersect(weight_gene, rownames(normalize_data))
    wg = gw[gs] # 包含权重和基因名称的向量
    
    ici_expr = normalize_data[gs, indd, drop=F]
    if(zscore){
      ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
    }
    for(xi in rownames(ici_expr)){
      ici_expr[xi,] = ici_expr[xi,] * wg[xi]
    }
    
    ici_means = colMeans(ici_expr)
    rocc = roc(ici_group, ici_means, levels=c("Response", "Non_response"), ci=TRUE)
    return(rocc)
    
  }) -> cl
  
  pval = sapply(cl, get_pval_roc)
  return(pval)
  
}) -> pl

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                         计算不同phase的加权得分
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plist = list()
for(l in 1:2){
  
  ## (1).  计算得分返回分组
  gw = gene_weight[,l]
  weight_gene = names(gw)[which(abs(gw) > weight_cut)]
  gw2 = gw[weight_gene]
  mat = ct@MURP$data_matrix[, names(gw2)]
  mat2 = apply(mat, 1, function(x) x*gw2)
  s = colMeans(mat2)
  c1 = which(ct@MURP$murp_cellinfo[,paste0("C0_",l)]==1)
  c2 = which(ct@MURP$murp_cellinfo[,paste0("C0_",l)]==2)
  dff = data.frame(score = s[c(c1,c2)], 
                  group = c(rep("Phase 1", length(c1)),
                            rep("Phase 2", length(c2))),
                  traj = paste0("Trajectory ",l))
  
  # dff = do.call(rbind, lx) %>% as.data.frame
  
  ## (2).  计算显著性
  annotation_df = get_pval_df(dff = dff, 
                              attr_name = c("traj"),
                              value_name = "score", 
                              group_name = "group", 
                              test_method = "t.test")
  annotation_df$label = sapply(annotation_df$pval, function(p){
    p2 = as.numeric(p)
    if(is.na(p2)){
      p
    }else{
      pval_lab(p2, lab = T)
    }
  })
  
  ## (3).  画图
  p <- ggplot(dff) + 
    geom_boxplot(aes(x = group, y = score, fill = group),
                 position = "dodge", width=0.4, size = 0.1,outlier.shape = NA) +
    labs(x = "", y = "Weighted mean of HWGs", fill = "") +
    # scale_fill_manual(values = c("#747474FF"))
    scale_fill_grey() +
    facet_wrap(~traj, nrow = 1) + 
    guides(fill = "none") + 
    # scale_y_continuous(limits = c(-0.1,0.1)) +
    geom_signif( data = annotation_df,
                 aes(xmin = X1, xmax = X2, annotations = TeX(label, output = "character"),
                     y_position = y_position-0.1),
                 tip_length = 0,
                 size = 0.06,
                 textsize = 5,
                 vjust = -0.1,
                 linetype = "dashed",
                 test = wilcox.test,
                 manual = TRUE,
                 parse = TRUE) +
    theme(panel.background = element_rect(fill='transparent', color="black", size = 0.5),
          strip.text = element_text(size = 13),
          strip.background = element_rect(colour = "black", fill = "transparent", size = 0.5),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.border = element_rect(fill='transparent', color='black', size = 0.5),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(vjust = -1.5, size = 13, colour = 'black'),
          axis.title.y = element_text(vjust = 1.5, size = 13, colour = 'black'),
          axis.text.y.left = element_text(vjust = 0, hjust = 1, size = 13, colour = 'black'),
          axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, size = 13, colour = 'black'),
          # axis.text.x = element_text(size = 3.7,colour = 'black'),
          strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")),
          strip.text.y = element_text(margin = margin(0,0.05,0,0.05, "cm")))
  plist = append(plist, list(p))
}
p = wrap_plots(plist, nrow = 1)
ggsave("phase_box.pdf", p, width = 6, height = 4) # gse99254

# ----

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                         高权重基因的加权轨迹得分
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# zscore = TRUE
# as <- c("Jung2019","Cho2020","Auslander2018", "Riaz2017","TCGA-SKCM")
# 
# lapply(seq(0.01,0.1,0.01), function(weight_cut){
#   
#   cat(weight_cut, "----------------\n")
#   
#   lapply(1:3, function(l){
#     gw = gene_weight[,l]
#     weight_gene = names(gw)[which(abs(gw) > weight_cut)]
#     
#     lapply(as, function(index){
#       
#       indd = which(clin$author==index)
#       ici_group = clin[indd, "Response"]
#       gs = intersect(weight_gene, rownames(normalize_data))
#       wg = gw[gs] # 包含权重和基因名称的向量
#       
#       ici_expr = normalize_data[gs, indd, drop=F]
#       if(zscore){
#         ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
#       }
#       for(xi in rownames(ici_expr)){
#         ici_expr[xi,] = ici_expr[xi,] * wg[xi]
#       }
#       
#       ici_means = colMeans(ici_expr)
#       rocc = roc(ici_group, ici_means, levels=c("Response", "Non_response"), ci=TRUE)
#       return(rocc)
#       
#     }) -> cl
#     
#     # lapply(as, function(index){
#     #   
#     #   indd = which(clin$author==index)
#     #   ici_group = clin[indd, "Response"]
#     #   gs = intersect(weight_gene, rownames(normalize_data))
#     #   ici_expr = normalize_data[gs, indd, drop=F]
#     #   if(zscore){
#     #     ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
#     #   }
#     #   ici_means = colMeans(ici_expr)
#     #   return(ici_means)
#     #   
#     # }) -> scorel
#     
#     names(cl) = as
#     glab = sapply(seq_along(cl), function(i){
#       paste0(names(cl)[i], " (AUC=",as.numeric(auc(cl[[i]])) %>% round(4),")")
#     }) %>% set_names(as)
#     
#     pROC::ggroc(cl, size = 1, aes=c("colour"), alpha = 1, legacy.axes = TRUE) +
#       geom_abline(slope=1, intercept=0, color = "#3e3e23ff", alpha = 0.8, size = 0.6, linetype = "dashed") +
#       labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", title = paste0("Trajectory ", l)) +
#       scale_color_aaas(name = "", breaks = names(glab), labels = glab) +
#       guides(linetype="none") +
#       theme(panel.background = element_rect(fill='transparent', color="black", size = 0.3),
#             strip.text = element_text(size = 10),
#             strip.background = element_rect(colour = "black", fill = "transparent", size = 0.25),
#             panel.grid.minor=element_blank(),
#             panel.grid.major=element_blank(),
#             panel.border = element_rect(fill='transparent', color='black', size = 0.2),
#             plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
#             plot.title = element_text(vjust = 0, size = 16, colour = 'black'),
#             legend.key = element_rect( fill = "white"),
#             legend.position = c(0.65,0.2),
#             legend.background = element_rect(fill='transparent'),
#             legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
#             legend.key.size = unit(0.6, "cm"),
#             axis.ticks = element_blank(),
#             axis.title.x = element_text(vjust = 0, size = 16, colour = 'black'),
#             axis.title.y = element_text(vjust = 2, size = 18, colour = 'black',margin = margin(r = 5)),
#             axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 16, colour = 'black', margin = margin(t = 6)),
#             axis.text.y = element_text(vjust = 0.5, size = 16, colour = 'black',margin = margin(r = 6)) ) -> p2
#     p2 = as.ggplot(as.grob(p2))
#     
#     return(p2)
#   }) -> pl
#   p = wrap_plots(pl, nrow = 1) +
#     plot_annotation(title = paste0("Weight CUT: ",weight_cut)) &
#     theme(plot.title = element_text(size = 14, face="bold"))
#   
#   return(p)
# }) -> plist
# 
# pdf(paste0("5_combine_plot/ici/ici_weight_gene_roc_zscore",zscore,"_trajectory_weightmeans.pdf"), width = 15, height = 5.5)
# plist
# dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#            weight gene 选择phase 1 -  phase 2 作为bifurcate score
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# zscore = TRUE
# as <- c("Jung2019","Cho2020","Auslander2018", "Riaz2017","TCGA-SKCM")
# # as = unique(clin$author)
# 
# lapply(seq(0.01,0.1,0.01), function(weight_cut){
#   
#   cat(weight_cut, "----------------\n")
#   
#   lapply(1:2, function(l){
#     cat(l, "\n")
#     gw = gene_weight[,l]
#     weight_gene = names(gw)[which(abs(gw) > weight_cut)]
#     data_matrix = ct@assay$data_matrix[,weight_gene]
#     
#     c = sdf[,paste0("C0_",l)]
#     c1 = which(c==1)
#     c2 = which(c==2)
#     phase1 = colMeans(data_matrix[c1, weight_gene])
#     phase2 = colMeans(data_matrix[c2, weight_gene])
#     dif = phase1 - phase2
#     p1 = phase1[which(dif>0)] %>% names ## 区分phase1 signature
#     p2 = phase2[which(dif<0)] %>% names## 区分phase2 signature
#     
#     lapply(as, function(index){
#       
#       indd = which(clin$author==index)
#       ici_group = clin[indd, "Response"]
#       gs1 = intersect(p1, rownames(normalize_data))
#       gs2 = intersect(p2, rownames(normalize_data))
#       ici_expr_1 = normalize_data[gs1, indd, drop=F]
#       ici_expr_2 = normalize_data[gs2, indd, drop=F]
#       if(zscore){
#         ici_expr_1 = apply(ici_expr_1, 1, function(x) (x-mean(x))/sd(x)); ici_expr_1 = t(ici_expr_1)
#         ici_expr_2 = apply(ici_expr_2, 1, function(x) (x-mean(x))/sd(x)); ici_expr_2 = t(ici_expr_2)
#       }
#       ici_means = colMeans(ici_expr_1) - colMeans(ici_expr_2)
#       rocc = roc(ici_group, ici_means, levels=c("Response", "Non_response"), ci=TRUE)
#       return(rocc)
#       
#     }) -> cl
#     names(cl) = as
#     glab = sapply(seq_along(cl), function(i){
#       paste0(names(cl)[i], " (AUC=",as.numeric(auc(cl[[i]])) %>% round(4),")")
#     }) %>% set_names(as)
#     
#     pROC::ggroc(cl, size = 1, aes=c("colour"), alpha = 1, legacy.axes = TRUE) +
#       geom_abline(slope=1, intercept=0, color = "#3e3e23ff", alpha = 0.8, size = 0.6, linetype = "dashed") +
#       labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", title = paste0("Trajectory ", l)) +
#       scale_color_npg(name = "", breaks = names(glab), labels = glab) +
#       guides(linetype="none") +
#       theme(panel.background = element_rect(fill='transparent', color="black", size = 0.3),
#             strip.text = element_text(size = 10),
#             strip.background = element_rect(colour = "black", fill = "transparent", size = 0.25),
#             panel.grid.minor=element_blank(),
#             panel.grid.major=element_blank(),
#             panel.border = element_rect(fill='transparent', color='black', size = 0.2),
#             plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
#             plot.title = element_text(vjust = 0, size = 16, colour = 'black'),
#             legend.key = element_rect( fill = "white"),
#             legend.position = c(0.65,0.2),
#             legend.background = element_rect(fill='transparent'),
#             legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
#             legend.key.size = unit(0.6, "cm"),
#             axis.ticks = element_blank(),
#             axis.title.x = element_text(vjust = 0, size = 16, colour = 'black'),
#             axis.title.y = element_text(vjust = 2, size = 18, colour = 'black',margin = margin(r = 5)),
#             axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 16, colour = 'black', margin = margin(t = 6)),
#             axis.text.y = element_text(vjust = 0.5, size = 16, colour = 'black',margin = margin(r = 6)) ) -> p2
#     p2 = as.ggplot(as.grob(p2))
#     
#     return(p2)
#   }) -> pl
#   p = wrap_plots(pl, nrow = 1) +
#     plot_annotation(title = paste0("Weight CUT: ",weight_cut)) &
#     theme(plot.title = element_text(size = 14, face="bold"))
#   
#   return(p)
# }) -> plist
# 
# pdf(paste0("5_combine_plot/ici/ici_weight_gene_roc_zscore",zscore,"_trajectory_bifscore.pdf"), width = 10, height = 5.5)
# plist
# dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#      weight gene 选择phase 1 -  phase 2 作为bifurcate score - 森林图
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# as <- c("Jung2019","Cho2020","Auslander2018", "Riaz2017","TCGA-SKCM")
# zscore = TRUE
# cancer_types = c("Jung2019" = "NSCLC", 
#                  "Cho2020" = "NSCLC",
#                  "Auslander2018" = "Melanoma",
#                  "Hugo2016" = "Melanoma",
#                  "Amato2020" = "Melanoma",
#                  "TCGA-SKCM" = "SKCM", 
#                  "Thibaudin2023" = "CRC",
#                  "Riaz2017" = "Melanoma")
# as = unique(clin$author)
# lapply(seq(0.01,0.09,0.01), function(weight_cut){
#   
#   cat(weight_cut, "----------------\n")
#   
#   lapply(1:2, function(l){
#     cat(l, "\n")
#     gw = gene_weight[,l]
#     weight_gene = names(gw)[which(abs(gw) > weight_cut)]
#     data_matrix = ct@assay$data_matrix[,weight_gene]
#     
#     c = sdf[,paste0("C0_",l)]
#     c1 = which(c==1)
#     c2 = which(c==2)
#     phase1 = colMeans(data_matrix[c1, weight_gene])
#     phase2 = colMeans(data_matrix[c2, weight_gene])
#     dif = phase1 - phase2
#     p1 = phase1[which(dif>0)] %>% names ## 区分phase1 signature
#     p2 = phase2[which(dif<0)] %>% names ## 区分phase2 signature
#     
#     lapply(as, function(index){
#       
#       indd = which(clin$author==index)
#       ici_group = clin[indd, "Response"]
#       gs1 = intersect(p1, rownames(normalize_data))
#       gs2 = intersect(p2, rownames(normalize_data))
#       ici_expr_1 = normalize_data[gs1, indd, drop=F]
#       ici_expr_2 = normalize_data[gs2, indd, drop=F]
#       if(zscore){
#         ici_expr_1 = apply(ici_expr_1, 1, function(x) (x-mean(x))/sd(x)); ici_expr_1 = t(ici_expr_1)
#         ici_expr_2 = apply(ici_expr_2, 1, function(x) (x-mean(x))/sd(x)); ici_expr_2 = t(ici_expr_2)
#       }
#       ici_means = colMeans(ici_expr_1) - colMeans(ici_expr_2)
#       
#       data = data.frame(group = factor(ici_group, levels = c("Response", "Non_response")), 
#                         ici_means = ici_means)
#       
#       
#       mod <- glm(formula = group ~ ici_means, data = data, family = "binomial")
#       mod_df <- cbind(exp(cbind(coef(mod), confint(mod))), 
#                       summary(mod)$coefficients[,4])
#       colnames(mod_df) <- c("OR", "2.5 %", "97.5 %", "P-value")
#       mod_df <- mod_df[2,,drop=F]
#       mod_df <- data.frame(mod_df, 
#                            lab = paste0( mod_df[,1] %>% round(3), "(", 
#                                          mod_df[,2] %>% round(3), "-",
#                                          mod_df[,3] %>% round(3), ")"),
#                            author = index, 
#                            traj = l, 
#                            weight_cut = weight_cut)
#       return(mod_df)
#       
#     })  %>% do.call(rbind, .)
#   })  %>% do.call(rbind, .) -> dff
#   
#   colnames(dff)[1:4] = c("mean", "lower", "upper", "pval")
#   dff$pval = ifelse(as.numeric(dff$pval) < 0.001, "<0.001", sprintf("%.3f", dff$pval))
#   dff$cancer = cancer_types[dff$author]
#   dff[,1:3] = log(dff[,1:3])
#   
#   ggplot(dff) +
#     geom_errorbar(aes(x = reorder(author, mean ), y = mean,
#                       ymin = lower, ymax = upper ),
#                   width = 0.4) +
#     geom_point(aes(x = reorder(author, desc(mean) ), y = mean),
#                color='darkorange3', size = 2.5) +
#     facet_wrap(~traj) +
#     coord_flip() +
#     geom_hline(yintercept = 0.5,linetype = "dashed") +
#     labs(x = "", y = "log(Odds Ratio)") +
#     theme(panel.background = element_rect(fill='transparent', color="black"),
#           strip.background = element_rect(fill='transparent', color="black"),
#           plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
#           plot.title = element_text(size = 15, hjust = 0),
#           panel.grid.minor=element_blank(),
#           panel.grid.major=element_blank(),
#           panel.border = element_rect(fill='transparent', color='black'),
#           legend.key = element_rect( fill = "white"),
#           legend.background = element_rect(fill='transparent'),
#           legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
#           legend.title = element_text(vjust = 0.4, size = 13, colour = 'black'),
#           legend.key.size = unit(0.6, "cm"),
#           legend.justification = "center",
#           axis.ticks = element_blank(),
#           axis.title.x = element_text(vjust = -1.5, size = 16, colour = 'black'),
#           axis.title.y = element_text(vjust = 1.5, size = 16, colour = 'black'),
#           axis.text.x = element_text(vjust = 1, size = 14,colour = 'black'),
#           axis.text.y = element_text(vjust = 1, size = 14,colour = 'black')) -> p
#   return(p)
#   
# }) -> pl
# 
# pdf(paste0("5_combine_plot/ici/ici_weight_gene_roc_zscore",zscore,"_trajectory_bifscore_pointrange.pdf"), width = 6, height = 3)
# pl
# dev.off()

# pdf(file=paste0("5_combine_plot/ici/1.pdf"), width=20, height=16, onefile = FALSE)
# dff %>%
#   forestplot::forestplot(labeltext = c(traj, lab, author),
#                          # clip = c(-2,80),
#                          clip = c(1,100),
#                          lty.ci = "solid",
#                          lwd.zero = 0.5,
#                          lwd.ci = 1.5,                 # 误差条的线的宽度
#                          ci.vertices.height = 0.1,     # 误差条末端的长度
#                          fn.ci_norm = fpDrawCircleCI,
#                          zero = 1,
#                          # xlog = T,
#                          align = c("c"),
#                          boxsize = 0.2,
#                          line.margin = 10,
#                          lines = unit(2,'cm'),
#                          lineheight = unit(10,'mm'),    # 设置图形中的行距
#                          colgap = unit(10,'mm'),        # 设置图形中的列间距
#                          title = "",
#                          xlab = "Odds ratio",
#                          txt_gp = fpTxtGp(ticks = gpar(cex=1),xlab = gpar(cex = 1.4)) ) %>%
#   fp_set_style(box = "red", summary = "#8B008B", lines = "black",zero = "black") %>%
#   fp_add_header(traj = c("","Trajectory"),
#                 # weight_cut = c("", "Gene Weight \nCutoff"),
#                 lab = c("", "Odds ratio\n95% CI"),
#                 pval = c(" ", "P value"),
#                 author = c(" ", "Author")) %>%
#   fp_set_zebra_style("grey90") %>% print
# dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                             合并3个数据的ici
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# weight_cut = 0.05
# indexs = c("Cho2020", "Auslander2018", "TCGA-SKCM")
# 
# lapply(1:3, function(l){
#   lapply(indexs, function(index){
#     indd = which(clin$author==index)
#     ici_group = clin[indd, "Response"]
#     ici_expr = normalize_data[, indd, drop=F]
#     if(zscore){
#       ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
#     }
#     gw = gene_weight[,l]
#     weight_gene = names(gw)[which(abs(gw) > weight_cut)]
#     weight_gene = intersect(weight_gene, rownames(ici_expr))
#     
#     ici_means = apply(ici_expr[weight_gene,,drop=F], 2, mean)
#     rocc = roc(ici_group, ici_means, levels=c("Response", "Non_response"), ci=TRUE) 
#     return(rocc)
#   }) -> xt
#   names(xt) <- indexs
#   return(xt)
# }) -> lx
# 
# plist = list()
# for(l in 1:3){
#   pROC::ggroc(list("Auslander2018" = lx[[l]][[2]], "Cho2020" = lx[[l]][[1]], "TCGA-SKCM" = lx[[l]][[3]]),
#               size = 1, aes=c("colour"), alpha = 1, legacy.axes = TRUE) +
#     geom_abline(slope=1, intercept=0, color = "#3e3e23ff", alpha = 0.8, size = 0.6, linetype = "dashed") +
#     labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") +
#     scale_color_aaas(name = "", breaks = c("Auslander2018", "Cho2020", "TCGA-SKCM"),
#                      labels = c(paste0("Auslander2018 (AUC=",as.numeric(auc( lx[[l]][[2]] )) %>% round(2),")"),
#                                 paste0("Cho2020 (AUC=",as.numeric(auc( lx[[l]][[1]] )) %>% round(2),")"),
#                                 paste0("TCGA-SKCM (AUC=",as.numeric(auc( lx[[l]][[3]] )) %>% round(2),")")  )) +
#     guides(linetype="none") +
#     theme(panel.background = element_rect(fill='transparent', color="black", size = 0.3),
#           strip.text = element_text(size = 10),
#           strip.background = element_rect(colour = "black", fill = "transparent", size = 0.25),
#           panel.grid.minor=element_blank(),
#           panel.grid.major=element_blank(),
#           panel.border = element_rect(fill='transparent', color='black', size = 0.2),
#           plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
#           plot.title = element_text(vjust = 0, size = 16, colour = 'black'),
#           legend.key = element_rect( fill = "white"),
#           legend.position = c(0.65,0.2),
#           legend.background = element_rect(fill='transparent'),
#           legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
#           legend.key.size = unit(0.6, "cm"),
#           axis.ticks = element_blank(),
#           axis.title.x = element_text(vjust = 0, size = 16, colour = 'black'),
#           axis.title.y = element_text(vjust = 2, size = 18, colour = 'black',margin = margin(r = 5)),
#           axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 16, colour = 'black', margin = margin(t = 6)),
#           axis.text.y = element_text(vjust = 0.5, size = 16, colour = 'black',margin = margin(r = 6)) ) -> p2
#   plist = append(plist, list(p2))
# }
# 
# pdf(paste0("5_combine_plot/ici/ici_auc_weight_", weight_cut, "_zscore_", zscore, "_combine.pdf"), width = 4.6, height = 4.5)
# plist
# dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                       合并3个数据的样本画一张图
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# weight_cut = 0.05
# indexs = c("Cho2020", "Auslander2018", "TCGA-SKCM")
# indd = which(clin$author %in% indexs)
# # indd = 1:nrow(clin)
# ici_group = clin[indd, "Response"]
# ici_expr = normalize_data[, indd, drop=F]
# if(zscore){
#   ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
# }
# lapply(1:3, function(l){
#   gw = gene_weight[,l]
#   weight_gene = names(gw)[which(abs(gw) > weight_cut)]
#   weight_gene = intersect(weight_gene, rownames(ici_expr))
#   
#   ici_means = apply(ici_expr[weight_gene,,drop=F], 2, mean)
#   rocc = roc(ici_group, ici_means, levels=c("Response", "Non_response"), ci=TRUE)
#   return(rocc)
# }) -> l
# 
# pROC::ggroc(list("Trajectory 1" = l[[1]], "Trajectory 2" = l[[2]], "Trajectory 3" = l[[3]]),
#             size = 1.5, aes=c("colour"), alpha = 1, legacy.axes = TRUE) +
#   geom_abline(slope=1, intercept=0, color = "#3e3e23ff", alpha = 0.8, size = 0.6, linetype = "dashed") +
#   labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") +
#   scale_color_aaas(name = "", breaks = c("Trajectory 1", "Trajectory 2", "Trajectory 3"),
#                    labels = c(paste0("Trajectory 1 (AUC=",as.numeric(auc(l[[1]])) %>% round(4),")"),
#                               paste0("Trajectory 2 (AUC=",as.numeric(auc(l[[2]])) %>% round(4),")"),
#                               paste0("Trajectory 3 (AUC=",as.numeric(auc(l[[3]])) %>% round(4),")")  )) +
#   guides(linetype="none") +
#   theme(panel.background = element_rect(fill='transparent', color="black", size = 0.3),
#         strip.text = element_text(size = 10),
#         strip.background = element_rect(colour = "black", fill = "transparent", size = 0.25),
#         panel.grid.minor=element_blank(),
#         panel.grid.major=element_blank(),
#         panel.border = element_rect(fill='transparent', color='black', size = 0.2),
#         plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
#         plot.title = element_text(vjust = 0, size = 16, colour = 'black'),
#         legend.key = element_rect( fill = "white"),
#         legend.position = c(0.65,0.2),
#         legend.background = element_rect(fill='transparent'),
#         legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
#         legend.key.size = unit(0.6, "cm"),
#         axis.ticks = element_blank(),
#         axis.title.x = element_text(vjust = 0, size = 16, colour = 'black'),
#         axis.title.y = element_text(vjust = 2, size = 18, colour = 'black',margin = margin(r = 5)),
#         axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 16, colour = 'black', margin = margin(t = 6)),
#         axis.text.y = element_text(vjust = 0.5, size = 16, colour = 'black',margin = margin(r = 6)) ) -> p2
# pdf(paste0("5_combine_plot/ici/ici_weight_gene_roc_zscore",zscore,"_allsample_3datasets.pdf"), width = 5, height = 4.7)
# p2
# dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                 AUC的柱状图
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# lapply(seq(0.01,0.1,0.01), function(weight_cut){
#   
#   cat(weight_cut, "----------------\n")
#   
#   lapply(unique(clin$author), function(index){
#     
#     indd = which(clin$author==index)
#     ici_group = clin[indd, "Response"]
#     ici_expr = normalize_data[, indd, drop=F]
#     if(zscore){
#       ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
#     }
#     
#     lapply(1:3, function(l){
#       gw = gene_weight[,l]
#       weight_gene = names(gw)[which(abs(gw) > weight_cut)]
#       weight_gene = intersect(weight_gene, rownames(ici_expr))
#       
#       ici_means = apply(ici_expr[weight_gene,,drop=F], 2, mean)
#       rocc = roc(ici_group, ici_means, levels=c("Response", "Non_response"), ci=TRUE) 
#       
#       df = data.frame(weight_cut = weight_cut,
#                       auc = as.numeric(auc(rocc)) %>% round(4),
#                       traj = l,
#                       author = index)
#       return(df)
#     })  %>% do.call(rbind, .)
#   }) %>% do.call(rbind, .)
# }) %>% do.call(rbind, .) -> aucl

## 柱状图
# weight_cutf = 0.05
# aucl = aucl %>% filter(weight_cut == weight_cutf, author %in% c("Cho2020", "Auslander2018", "TCGA-SKCM"))
# aucl$traj = factor(aucl$traj, levels = 1:3)
# ggplot(aucl, aes(x = author, y = auc, group = traj)) + 
#   geom_col(aes(fill = traj), width = 0.8, position = "dodge") +
#   # geom_text(aes(label = pval_lab2),  parse = T, 
#   #           position = position_dodge(width = .9), vjust = 0.2, size = 10 / .pt) +
#   labs(x = "", y = "AUC", fill = "") +
#   # facet_wrap(~ weight_cut, scale = "free") +
#   scale_y_continuous(limits = c(0,1)) +
#   scale_fill_grey() +
#   theme(panel.background = element_rect(fill='transparent', color="black", size = 0.2),
#         strip.text = element_text(size = 10),
#         strip.background = element_rect(colour = "black", fill = "transparent", size = 0.25),
#         panel.grid.minor=element_blank(),
#         panel.grid.major=element_blank(),
#         panel.border = element_rect(fill='transparent', color='black', size = 0.2),
#         plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
#         axis.ticks = element_blank(),
#         axis.title.x = element_text(vjust = -1.5, size = 10, colour = 'black'),
#         axis.title.y = element_text(vjust = 1.5, size = 10, colour = 'black'),
#         axis.text.y.left = element_text(vjust = 0, hjust = 1, size = 10, colour = 'black'),
#         axis.text.x = element_text(vjust = 1.03, hjust = 1, angle = 45, size = 10, colour = 'black'),
#         # axis.text.x = element_text(size = 4,colour = 'black'),
#         strip.text.x = element_text(margin = margin(0.05, 0, 0.05 ,0, "cm")),
#         strip.text.y = element_text(margin = margin(0, 0.05, 0, 0.05, "cm"))) -> p
# ggsave(paste0("5_combine_plot/ici/ici_auc_weight_", weight_cutf, "_zscore_", zscore, "_bar.pdf"), p, width = 3.5, height = 3.3)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                        导入ICI数据 1个数据包含3个轨迹
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# as <- c("Jung2019","Cho2020","Auslander2018", "Riaz2017","Hugo2016","TCGA-SKCM")
# 
# lapply(seq(0.01,0.1,0.01), function(weight_cut){
#   
#   cat(weight_cut, "----------------\n")
#   
#   lapply(as, function(index){
#     
#     indd = which(clin$author==index)
#     ici_group = clin[indd, "Response"]
#     ici_expr = normalize_data[, indd, drop=F]
#     if(zscore){
#       ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
#     }
#     
#     lapply(1:3, function(l){
#       gw = gene_weight[,l]
#       weight_gene = names(gw)[which(abs(gw) > weight_cut)]
#       weight_gene = intersect(weight_gene, rownames(ici_expr))
#       
#       ici_means = apply(ici_expr[weight_gene,,drop=F], 2, mean)
#       rocc = roc(ici_group, ici_means, levels=c("Response", "Non_response"), ci=TRUE)
#       return(rocc)
#     }) -> l
#     pROC::ggroc(list("Trajectory 1" = l[[1]], "Trajectory 2" = l[[2]], "Trajectory 3" = l[[3]]),
#                 size = 1.5, aes=c("colour"), alpha = 1, legacy.axes = TRUE) +
#       geom_abline(slope=1, intercept=0, color = "#3e3e23ff", alpha = 0.8, size = 0.6, linetype = "dashed") +
#       labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", title = paste0(index)) +
#       scale_color_aaas(name = "", breaks = c("Trajectory 1", "Trajectory 2", "Trajectory 3"),
#                        labels = c(paste0("Trajectory 1 (AUC=",as.numeric(auc(l[[1]])) %>% round(4),")"),
#                                   paste0("Trajectory 2 (AUC=",as.numeric(auc(l[[2]])) %>% round(4),")"),
#                                   paste0("Trajectory 3 (AUC=",as.numeric(auc(l[[3]])) %>% round(4),")")  )) +
#       guides(linetype="none") +
#       theme(panel.background = element_rect(fill='transparent', color="black", size = 0.3),
#             strip.text = element_text(size = 10),
#             strip.background = element_rect(colour = "black", fill = "transparent", size = 0.25),
#             panel.grid.minor=element_blank(),
#             panel.grid.major=element_blank(),
#             panel.border = element_rect(fill='transparent', color='black', size = 0.2),
#             plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
#             plot.title = element_text(vjust = 0, size = 16, colour = 'black'),
#             legend.key = element_rect( fill = "white"),
#             legend.position = c(0.65,0.2),
#             legend.background = element_rect(fill='transparent'),
#             legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
#             legend.key.size = unit(0.6, "cm"),
#             axis.ticks = element_blank(),
#             axis.title.x = element_text(vjust = 0, size = 16, colour = 'black'),
#             axis.title.y = element_text(vjust = 2, size = 18, colour = 'black',margin = margin(r = 5)),
#             axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 16, colour = 'black', margin = margin(t = 6)),
#             axis.text.y = element_text(vjust = 0.5, size = 16, colour = 'black',margin = margin(r = 6)) ) -> p2
#     p2 = as.ggplot(as.grob(p2))
#     return(p2)
#   }) -> plist
#   p = wrap_plots(plist, nrow = 1) +
#     plot_annotation(title = paste0("Weight CUT: ",weight_cut)) &
#     theme(plot.title = element_text(size = 14, face="bold"))
#   
#   return(p)
# }) -> plist
# 
# pdf(paste0("5_combine_plot/ici/ici_weight_gene_roc_zscore",zscore,"_2.pdf"), width = 41, height = 5.5)
# plist
# dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                           1个轨迹包含多个数据 ***
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# zscore = TRUE
# as <- c("Jung2019","Cho2020","Auslander2018", "Riaz2017","TCGA-SKCM")
# 
# lapply(seq(0.01,0.1,0.01), function(weight_cut){
#   
#   cat(weight_cut, "----------------\n")
#   
#   lapply(1:3, function(l){
#     gw = gene_weight[,l]
#     weight_gene = names(gw)[which(abs(gw) > weight_cut)]
#     
#     lapply(as, function(index){
#       
#       indd = which(clin$author==index)
#       ici_group = clin[indd, "Response"]
#       gs = intersect(weight_gene, rownames(normalize_data))
#       ici_expr = normalize_data[gs, indd, drop=F]
#       if(zscore){
#         ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
#       }
#       ici_means = colMeans(ici_expr)
#       rocc = roc(ici_group, ici_means, levels=c("Response", "Non_response"), ci=TRUE)
#       return(rocc)
#       
#     }) -> cl
#     
#     lapply(as, function(index){
#       
#       indd = which(clin$author==index)
#       ici_group = clin[indd, "Response"]
#       gs = intersect(weight_gene, rownames(normalize_data))
#       ici_expr = normalize_data[gs, indd, drop=F]
#       if(zscore){
#         ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
#       }
#       ici_means = colMeans(ici_expr)
#       return(ici_means)
#     }) -> scorel
#     
#     names(cl) = as
#     glab = sapply(seq_along(cl), function(i){
#       paste0(names(cl)[i], " (AUC=",as.numeric(auc(cl[[i]])) %>% round(4),")")
#     }) %>% set_names(as)
#     
#     pROC::ggroc(cl, size = 1, aes=c("colour"), alpha = 1, legacy.axes = TRUE) +
#       geom_abline(slope=1, intercept=0, color = "#3e3e23ff", alpha = 0.8, size = 0.6, linetype = "dashed") +
#       labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", title = paste0("Trajectory ", l)) +
#       scale_color_aaas(name = "", breaks = names(glab), labels = glab) +
#       guides(linetype="none") +
#       theme(panel.background = element_rect(fill='transparent', color="black", size = 0.3),
#             strip.text = element_text(size = 10),
#             strip.background = element_rect(colour = "black", fill = "transparent", size = 0.25),
#             panel.grid.minor=element_blank(),
#             panel.grid.major=element_blank(),
#             panel.border = element_rect(fill='transparent', color='black', size = 0.2),
#             plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
#             plot.title = element_text(vjust = 0, size = 16, colour = 'black'),
#             legend.key = element_rect( fill = "white"),
#             legend.position = c(0.65,0.2),
#             legend.background = element_rect(fill='transparent'),
#             legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
#             legend.key.size = unit(0.6, "cm"),
#             axis.ticks = element_blank(),
#             axis.title.x = element_text(vjust = 0, size = 16, colour = 'black'),
#             axis.title.y = element_text(vjust = 2, size = 18, colour = 'black',margin = margin(r = 5)),
#             axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 16, colour = 'black', margin = margin(t = 6)),
#             axis.text.y = element_text(vjust = 0.5, size = 16, colour = 'black',margin = margin(r = 6)) ) -> p2
#     p2 = as.ggplot(as.grob(p2))
#     
#     return(p2)
#   }) -> pl
#   p = wrap_plots(pl, nrow = 1) +
#     plot_annotation(title = paste0("Weight CUT: ",weight_cut)) &
#     theme(plot.title = element_text(size = 14, face="bold"))
#   
#   return(p)
# }) -> plist
# 
# pdf(paste0("5_combine_plot/ici/ici_weight_gene_roc_zscore",zscore,"_trajectory.pdf"), width = 15, height = 5.5)
# plist
# dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                           1个轨迹包含多个数据 ***
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# zscore = TRUE
# as <- c("Jung2019","Cho2020","Auslander2018", "Riaz2017","TCGA-SKCM")
# 
# lapply(seq(0.01,0.1,0.01), function(weight_cut){
#   
#   cat(weight_cut, "----------------\n")
#   
#   lapply(1:3, function(l){
#     gw = gene_weight[,l]
#     weight_gene = names(gw)[which(abs(gw) > weight_cut)]
#     
#     lapply(as, function(index){
#       
#       indd = which(clin$author==index)
#       ici_group = clin[indd, "Response"]
#       gs = intersect(weight_gene, rownames(normalize_data))
#       ici_expr = normalize_data[gs, indd, drop=F]
#       if(zscore){
#         ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
#       }
#       ici_means = colMeans(ici_expr)
#       rocc = roc(ici_group, ici_means, levels=c("Response", "Non_response"), ci=TRUE)
#       return(rocc)
#       
#     }) -> cl
#     names(cl) = as
#     glab = sapply(seq_along(cl), function(i){
#       paste0(names(cl)[i], " (AUC=",as.numeric(auc(cl[[i]])) %>% round(4),")")
#     }) %>% set_names(as)
#     
#     pROC::ggroc(cl, size = 1, aes=c("colour"), alpha = 1, legacy.axes = TRUE) +
#       geom_abline(slope=1, intercept=0, color = "#3e3e23ff", alpha = 0.8, size = 0.6, linetype = "dashed") +
#       labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", title = paste0("Trajectory ", l)) +
#       scale_color_npg(name = "", breaks = names(glab), labels = glab) +
#       guides(linetype="none") +
#       theme(panel.background = element_rect(fill='transparent', color="black", size = 0.3),
#             strip.text = element_text(size = 10),
#             strip.background = element_rect(colour = "black", fill = "transparent", size = 0.25),
#             panel.grid.minor=element_blank(),
#             panel.grid.major=element_blank(),
#             panel.border = element_rect(fill='transparent', color='black', size = 0.2),
#             plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
#             plot.title = element_text(vjust = 0, size = 16, colour = 'black'),
#             legend.key = element_rect( fill = "white"),
#             legend.position = c(0.65,0.2),
#             legend.background = element_rect(fill='transparent'),
#             legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
#             legend.key.size = unit(0.6, "cm"),
#             axis.ticks = element_blank(),
#             axis.title.x = element_text(vjust = 0, size = 16, colour = 'black'),
#             axis.title.y = element_text(vjust = 2, size = 18, colour = 'black',margin = margin(r = 5)),
#             axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 16, colour = 'black', margin = margin(t = 6)),
#             axis.text.y = element_text(vjust = 0.5, size = 16, colour = 'black',margin = margin(r = 6)) ) -> p2
#     p2 = as.ggplot(as.grob(p2))
#     
#     return(p2)
#   }) -> pl
#   p = wrap_plots(pl, nrow = 1) +
#     plot_annotation(title = paste0("Weight CUT: ",weight_cut)) &
#     theme(plot.title = element_text(size = 14, face="bold"))
#   
#   return(p)
# }) -> plist
# 
# pdf(paste0("5_combine_plot/ici/ici_weight_gene_roc_zscore",zscore,"_trajectory.pdf"), width = 15, height = 5.5)
# plist
# dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                     weight gene 选择phase 1上调的
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# zscore = TRUE
# as <- c("Jung2019","Cho2020","Auslander2018", "Riaz2017","TCGA-SKCM")
# # as = unique(clin$author)
# lapply(seq(0.01,0.1,0.01), function(weight_cut){
#   
#   cat(weight_cut, "----------------\n")
#   
#   lapply(1:2, function(l){
#     cat(l, "\n")
#     gw = gene_weight[,l]
#     weight_gene = names(gw)[which(abs(gw) > weight_cut)]
#     data_matrix = ct@assay$data_matrix[,weight_gene]
#     
#     c = sdf[,paste0("C0_",l)]
#     c1 = which(c==1)
#     c2 = which(c==2)
#     phase1 = colMeans(data_matrix[c1, weight_gene])
#     phase2 = colMeans(data_matrix[c2, weight_gene])
#     dif = phase1 - phase2
#     p1 = phase1[which(dif>0)] %>% names ## 区分phase1 signature
#     p2 = phase2[which(dif<0)] %>% names## 区分phase2 signature
#     
#     lapply(as, function(index){
#       
#       indd = which(clin$author==index)
#       ici_group = clin[indd, "Response"]
#       gs = intersect(p1, rownames(normalize_data))
#       ici_expr = normalize_data[gs, indd, drop=F]
#       if(zscore){
#         ici_expr = apply(ici_expr, 1, function(x) (x-mean(x))/sd(x)); ici_expr = t(ici_expr)
#       }
#       ici_means = colMeans(ici_expr)
#       rocc = roc(ici_group, ici_means, levels=c("Response", "Non_response"), ci=TRUE)
#       return(rocc)
#       
#     }) -> cl
#     names(cl) = as
#     glab = sapply(seq_along(cl), function(i){
#       paste0(names(cl)[i], " (AUC=",as.numeric(auc(cl[[i]])) %>% round(4),")")
#     }) %>% set_names(as)
#     
#     pROC::ggroc(cl, size = 1, aes=c("colour"), alpha = 1, legacy.axes = TRUE) +
#       geom_abline(slope=1, intercept=0, color = "#3e3e23ff", alpha = 0.8, size = 0.6, linetype = "dashed") +
#       labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", title = paste0("Trajectory ", l)) +
#       scale_color_npg(name = "", breaks = names(glab), labels = glab) +
#       guides(linetype="none") +
#       theme(panel.background = element_rect(fill='transparent', color="black", size = 0.3),
#             strip.text = element_text(size = 10),
#             strip.background = element_rect(colour = "black", fill = "transparent", size = 0.25),
#             panel.grid.minor=element_blank(),
#             panel.grid.major=element_blank(),
#             panel.border = element_rect(fill='transparent', color='black', size = 0.2),
#             plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
#             plot.title = element_text(vjust = 0, size = 16, colour = 'black'),
#             legend.key = element_rect( fill = "white"),
#             legend.position = c(0.65,0.2),
#             legend.background = element_rect(fill='transparent'),
#             legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
#             legend.key.size = unit(0.6, "cm"),
#             axis.ticks = element_blank(),
#             axis.title.x = element_text(vjust = 0, size = 16, colour = 'black'),
#             axis.title.y = element_text(vjust = 2, size = 18, colour = 'black',margin = margin(r = 5)),
#             axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 16, colour = 'black', margin = margin(t = 6)),
#             axis.text.y = element_text(vjust = 0.5, size = 16, colour = 'black',margin = margin(r = 6)) ) -> p2
#     p2 = as.ggplot(as.grob(p2))
#     
#     return(p2)
#   }) -> pl
#   p = wrap_plots(pl, nrow = 1) +
#     plot_annotation(title = paste0("Weight CUT: ",weight_cut)) &
#     theme(plot.title = element_text(size = 14, face="bold"))
#   
#   return(p)
# }) -> plist
# 
# pdf(paste0("5_combine_plot/ici/ici_weight_gene_roc_zscore",zscore,"_trajectory_p1.pdf"), width = 10, height = 5.5)
# plist
# dev.off()

