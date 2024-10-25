
julia_home = "/public/home/renjun/tool/julia-1.6.6/bin"
celltrek_root <- "/share/data6/tmp/renjun/CellTrekResult/CellTrek/celltrekr/R"
invisible(sapply(list.files(celltrek_root, full.names = T), source))
celltrek_root <- "/share/data6/tmp/renjun/CellTrekResult/CellTrek/celltrek_code_all"
source(paste0(celltrek_root,"/prepare.R"))

library(GeneOverlap)
library(purrr)

dir.create("5_combine_plot/ici")
setwd("/share/data6/tmp/renjun/CellTrekResult/CellTrek/celltrek_code_all/ici/")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                         ctrDB - 读取所有数据
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## count data
# dataset = c("Jung2019", "Cho2020", "Auslander2018", "Riaz2017", "Hugo2016","Amato2020", "TCGA-SKCM", "Thibaudin2023")

## normalize data
# dataset = c("Cui2021", "Amato2020", "Gide2019", "Hugo2016", "Jung2019", "Liu2019", "Riaz2017", "Thibaudin2023")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                   limma 
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

adjust_method = "BH"
dirn = paste0("limma_", adjust_method)
dir.create(dirn)

# 分数据集定义ICI
for(ds in dataset){
  
  expr = read.csv(paste0("/share/data6/tmp/renjun/Reference/ICI_RNA/summary_data/", ds, "/normalize_data.csv"), row.names = 1, header = T)
  clin = read.csv(paste0("/share/data6/tmp/renjun/Reference/ICI_RNA/summary_data/", ds, "/clin.csv"), row.names = 1, header = T)
  group = factor(clin$Response, c("Response", "Non_response"))
  
  library(limma)
  design = model.matrix(~ 0 + group)
  colnames(design) = levels(factor(group))
  rownames(design) = colnames(expr)
  
  # v = voom(expr)
  # v = voom(expr, design)
  # expr = log(expr+1)
  contrast = makeContrasts("Response-Non_response",levels=design)
  vfit = lmFit(expr, design)
  vfit = contrasts.fit(vfit, contrasts = contrast)
  efit = eBayes(vfit)
  deg = topTable(efit, number=200000, adjust.method="BH")
  colnames(deg)[c(1,4,5)] = c("logFC", "pval", "p.adj")
  
  x = deg %>% dplyr::filter(p.adj < 0.05, abs(logFC)>2) %>% dim
  cat(x, "\n")
  write.csv(deg, file = paste0(dirn,"/",ds,".csv"))
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                  deseq2
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

adjust_method = "BH"
dirn = paste0("deseq2_", adjust_method)
dir.create(dirn)

##  分数据集定义ICI
for(ds in dataset){
 
  cat(ds, "\n")
  
  rm(expr, clin, gene, deg, normalize_data, dds)
  expr = read.csv(paste0("/share/data6/tmp/renjun/Reference/ICI_RNA/summary_data/", ds, "/count.csv"), row.names = 1, header = T)
  clin = read.csv(paste0("/share/data6/tmp/renjun/Reference/ICI_RNA/summary_data/", ds, "/clin.csv"), row.names = 1, header = T)
  
  gene = rownames(expr)[which(apply(expr, 1, sum)>1)]
  expr = expr[gene, ]
  group = factor(clin$Response, c("Response", "Non_response"))
  
  library(DESeq2)
  dds = DESeqDataSetFromMatrix(count = expr, colData = clin, design= ~ Response)
  dds = DESeq(dds)
  
  normalize_data = counts(dds, normalize = T)
  deg = results(dds, contrast = c("Response","Response","Non_response"), pAdjustMethod = adjust_method)
  deg = as.data.frame(deg)
  colnames(deg)[c(2,5,6)] = c("logFC", "pval", "p.adj")
  
  x = deg %>% dplyr::filter(p.adj < 0.05, abs(logFC)>2) %>% dim
  cat(x, "\n")
  write.csv(deg, file = paste0(dirn,"/",ds,".csv"))
  save(dds, normalize_data, file = paste0(dirn, "/", ds,".rda"))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                               wilcox test
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method = "wilcox.test"
adjust_method = "BH"
dirn = paste0(method, "_", adjust_method)
dir.create(dirn)

# 分数据集定义ICI
for(ds in dataset){
  cat(ds, "\n")
  
  rm(expr, clin, gene, deg)
  expr = read.csv(paste0("/share/data6/tmp/renjun/Reference/ICI_RNA/summary_data/", ds, "/normalize_data.csv"), row.names = 1, header = T)
  clin = read.csv(paste0("/share/data6/tmp/renjun/Reference/ICI_RNA/summary_data/", ds, "/clin.csv"), row.names = 1, header = T)
  
  gene = rownames(expr)
  group = factor(clin$Response, c("Response", "Non_response"))
  
  lapply(gene, function(g){
    # cat(g, "\n")
    df = data.frame(group = group,
                    g = unlist(expr[g,]))
    pval = tryCatch(
      {
        if(method=="wilcox.test"){
          mod = wilcox.test(df$g~df$group)
          p = mod$p.value
        }
        if(method=="t.test"){
          mod = t.test(df$g~df$group)
          p = mod$p.value
        }
        p
      },
      error = function(e){
        NA
      }
    )
    c(gene = g, pval = pval)
  }) %>% do.call(rbind, .) %>% data.frame -> deg
  
  rownames(deg) = deg$gene
  deg[,2] = as.numeric(deg[,2])
  deg$p.adj = p.adjust(deg[,2], method = adjust_method)
  tapply(1:ncol(expr), group,
         function(x, y) apply(expr[,x,drop=F], 1, mean)) %>%
    do.call(cbind, .) -> g_means
  deg = cbind(deg, g_means)
  deg$logFC = log2(deg$Response / deg$Non_response)
  x = deg %>% dplyr::filter(pval < 0.05, abs(logFC)>2) %>% dim
  cat(dim(x), "\n")
  write.csv(deg, file = paste0(dirn,"/",ds,".csv"))
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                              fisher_test
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirn = "deseq2_BH"
# dirn = "wilcox.test_BH"
pval_cut = 0.05
fc = 2
dataset = list.files(paste0(celltrek_root, "/ici/", dirn))
dataset = gsub("\\.csv|\\.rda", "", dataset)
dataset = unique(dataset)
lapply(dataset, 
       function(index){
         ici = read.csv(paste0(celltrek_root, "/ici/", dirn, "/", index, ".csv"), row.names = 1)
         cat(dim(ici), "\n")
         # ici %>% dplyr::filter(pvalue < pval_cut, abs(logFC) > fc) %>% rownames
         ici %>% dplyr::filter(p.adj < pval_cut, abs(logFC) > fc) %>% rownames
       }) -> g_list
names(g_list) = dataset
source(paste0(celltrek_root, "/ici/ici_fisher_test.R"))