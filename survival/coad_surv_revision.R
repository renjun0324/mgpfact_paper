
# source("/share/data6/tmp/renjun/CellTrekResult/CellTrek/celltrek_code_all/survival/coad_surv_revision.R")

library(survival)
library(survminer)
library(philentropy)
library(forestplot)

xtcga = data.table::fread("./5_combine_plot/survival/TCGA.COAD.sampleMap_HiSeqV2_PANCAN.gz", data.table = F)
rownames(xtcga) = xtcga[,1]
xtcga = xtcga[,-1]
xclin = read.table("./5_combine_plot/survival/survival_COAD_survival.txt", sep = "\t", header = T)
rownames(xclin) = xclin$sample

x = intersect(colnames(xtcga), rownames(xclin))
tcga = xtcga[, x]
clin = xclin[x,]

# 定义结局
clin$status = clin$OS
clin$futime = clin$OS.time %>% as.numeric
clin$futime = clin$futime/365

# 去掉不包含生存时间的样本
tcga = tcga[,!is.na(clin$futime)]
clin = clin[!is.na(clin$futime),]
all.equal(rownames(clin),colnames(tcga))

# # 去掉生存时间为0的样本
tcga = tcga[,which(clin$futime!=0)]
clin = clin[which(clin$futime!=0),]
all.equal(rownames(clin),colnames(tcga))

# rmdx = c("TCGA-DM-A1D0-01",
#          "TCGA-DM-A1HA-01",
#          "TCGA-G4-6299-01",
#          "TCGA-DM-A0X9-01",
#          "TCGA-AZ-6601-01",
#          "TCGA-AD-6889-01",
#          "TCGA-DM-A28E-01",
#          "TCGA-G4-6304-01",
#          "TCGA-D5-5540-01",
#          "TCGA-D5-5541-01",
#          "TCGA-G4-6307-01",
#          "TCGA-A6-A56B-01",
#          "TCGA-G4-6626-01",
#          "TCGA-F4-6703-01",
#          "TCGA-F4-6855-01")
# idx = which(clin$sample %nin% rmdx)
# tcga = tcga[,idx]
# clin_tumor = clin[idx,]
clin_tumor = clin

sdf = ct@MetaData
gene_weight = GetGeneWeight(ct)

load(paste0("/share/data6/tmp/renjun/ExtraJobs/viper/1_tcga_rna_data/PanCanAtlas/COAD/clinical.rda"))
clin_tumor$tumor_stage = clinical[clin_tumor[,2],"ajcc_pathologic_tumor_stage"]
clin_tumor$age = clinical[clin_tumor[,2],"age_at_initial_pathologic_diagnosis"]
clin_tumor$gender = clinical[clin_tumor[,2],"gender"]
clin_tumor$purity = clinical[clin_tumor[,2],"purity"]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                               分轨迹计算
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plist <- list()
for(l in 1:2){
  
  ## (0). 保留权重基因
  
  gw = sort(abs(gene_weight[,l]), decreasing = T)
  genes = gw[1:ceiling(length(gw)*1)] %>% names
  # genes = gw[which(gw>cutf)] %>% names
  weight_gene = intersect(genes, rownames(tcga))
  
  tcga2 = tcga[weight_gene,,drop=F]
  data_matrix = ct@assay$data_matrix[,weight_gene]
  
  ## (1). 不区分up
  
  c = sdf[,paste0("C0_",l)]
  c1 = which(c==1)
  c2 = which(c==2)
  phase1 = colMeans(data_matrix[c1, weight_gene])
  phase2 = colMeans(data_matrix[c2, weight_gene])
  sapply(1:ncol(tcga2), function(i){
    sapply(list(phase1, phase2), function(x){
      cor(x, tcga2[names(x),i], method = "pearson")
    }) %>% which.max
  }) -> s_group
  
  ## (1). 利用phase1 up 和phase2 up
  
  # c = sdf[,paste0("C0_",l)]
  # c1 = which(c==1)
  # c2 = which(c==2)
  # phase1 = colMeans(data_matrix[c1, weight_gene])
  # phase2 = colMeans(data_matrix[c2, weight_gene])
  # dif = phase1 - phase2
  # p1 = phase1[which(dif>0)] ## 区分phase1 signature
  # p2 = phase2[which(dif<0)] ## 区分phase2 signature
  # sapply(1:ncol(tcga2), function(i){
  #   sapply(list(p1, p2), function(x){
  #     cor(x, tcga2[names(x),i], method = "pearson")
  #   }) %>% which.max
  # }) -> s_group
  
  if(dim(table(s_group))==1){
    p = ggplot()
    plist = append(plist, list(p))
    next
  }
  
  ## (4). 生存分析
  
  df = data.frame(clin_tumor, s_group = factor(paste0("g", s_group)) )
  fit <- survfit(Surv(futime, status) ~ s_group, data = df)
  
  ## (5). COX分析
  
  # cox = coxph(Surv(futime, status) ~ s_group, data = df)
  ## revision: 矫正掉混杂因素
  cox = coxph(Surv(futime, status) ~ s_group + age + gender + tumor_stage + purity, data = df)  

  coxSummary=summary(cox)
  coxP=coxSummary$coefficients[1,"Pr(>|z|)"]
  hr = round(coxSummary$conf.int[1,"exp(coef)"],3)
  lowerr = round(coxSummary$conf.int[1,"lower .95"], 3)
  highh = round(coxSummary$conf.int[1,"upper .95"],3)
  tx = paste0("HR = ", hr,
              "(", lowerr,
              "-", highh, ")",
              "\nlog rank P value = ", round(coxP,3),
              "\n")
  cat(tx,"\n")
  
  # (7). KM曲线
  legend.labs = c(paste0("Phase 1 (n=",fit$n[1],")"),
                  paste0("Phase 2 (n=", fit$n[2],")"))
  ggsurvplot(fit,
             data = df,
             font.main = c(14, "bold"),
             conf.int = FALSE,
             conf.int.style = "ribbon",
             conf.int.alpha = 0.2,
             font.x = 14,
             font.y = 14,
             font.tickslab = 12,
             # pval = TRUE,
             pval = tx,
             pval.size=4.5,
             # pval.coord=c(148, 0.25),
             pval.coord=c(1, 0.25),
             surv.median.line = "hv",
             size = 1.2,
             linetype = "strata",
             palette = "Dark2",
             legend = c(0.6, 0.83),
             legend.title = "",
             legend.labs = legend.labs,
             risk.table = F,
             ggtheme = rj.ftheme +
               theme(panel.grid = element_blank(),
                     legend.text = element_text(size = 12),
                     legend.title = element_text(size = 14)) ) +
    labs(x = "Overall Time (Years)",
         y = "Survival probability",
         title = paste0("TCGA tumors Trajectory ", l))-> p
  plist = append(plist, list(p$plot))
  
}

p = wrap_plots(plist, ncol = 3)
ggsave(paste0("coad_trajectory_survival_revision.pdf"), p, width = 14, height = 5, limitsize = FALSE)
