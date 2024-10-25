
setwd("/share/data6/tmp/renjun/CellTrekResult/CellTrekData/GSE123025/3_result/paper_for_monocle/5_combine_plot/originial_trajectory/data_diff_test_p7_ident_qval0.01_350**/")

BEAM_res = read.csv("BEAM_res_data.csv")
BEAM_res = BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
# beam.deg <- BEAM_res$gene_short_name[which(BEAM_res$qval<0.01)]
beam.deg <- BEAM_res$gene_short_name[1:500]

intersect(beam.deg, colnames(ct@assay$data_matrix))

intersect(beam.deg, weight_gene_list[[1]]) %>% length
intersect(beam.deg, weight_gene_list[[2]]) %>% length
intersect(beam.deg, weight_gene_list[[3]]) %>% length

intersect(intersect(beam.deg, colnames(ct@assay$data_matrix)), weight_gene_list[[1]]) %>% length
intersect(intersect(beam.deg, colnames(ct@assay$data_matrix)), weight_gene_list[[2]]) %>% length
intersect(intersect(beam.deg, colnames(ct@assay$data_matrix)), weight_gene_list[[3]]) %>% length

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                               
#                           和MGPfact的韦恩图
#                                                                              
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# gl <- append(weight_gene_list, list(beam.deg))
# ggVennDiagram(gl, label_alpha=0,label_size =3,
#               edge_size = 0.5,label ="count") +
#   scale_color_lancet()+
#   scale_fill_gradient(low="gray100",high = "gray95",guide="none")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                              
#                         准备转化之后的genelist
#                                                                              
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 1). get genes
genes = beam.deg
# genes = intersect(beam.deg, colnames(ct@assay$data_matrix))

## 2). bitr - gene_df 
tmp = bitr(genes, 
           fromType = "SYMBOL", 
           toType = c("ENSEMBL", "ENTREZID"),
           OrgDb = ifelse(species=="mouse", "org.Mm.eg.db", "org.Hs.eg.db"))
tmp <- dplyr::distinct(tmp,SYMBOL,.keep_all=TRUE) # 去重
colnames(tmp)[1] = "gene"
entrezid = tmp$ENTREZID

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                              
#                                 enrichment
#                                                                              
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 1. KEGG
x = enrichKEGG(entrezid,
               organism = ifelse(species=="mouse", "mmu", "hsa"),
               keyType = 'kegg', 
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               qvalueCutoff = 0.1,
               minGSSize = 10,
               maxGSSize = 300,
               use_internal_data = FALSE)
x = setReadable(x, OrgDb = ifelse(species=="mouse", "org.Mm.eg.db", "org.Hs.eg.db"), 
                keyType="ENTREZID") %>% data.frame
y = x %>% data.frame
# y = x[which(x$Count>10),] %>% data.frame
gr = stringr::str_split_fixed(y$GeneRatio,"/",2)
br = stringr::str_split_fixed(y$BgRatio,"/",2)
fe = (as.numeric(gr[,1])/as.numeric(gr[,2])) / (as.numeric(br[,1])/as.numeric(br[,2]))
y$fold_enrichment = fe
y$generatio = as.numeric(gr[,1])/as.numeric(gr[,2])
y$bgratio = as.numeric(br[,1])/as.numeric(br[,2])
kegg_result = y
save(kegg_result, file = paste0("enrich_kegg_list.rda"))

tryCatch({
    write.xlsx(kegg_result,
               file = paste0("enrich_kegg_combine.xlsx"),
               sheetName = "kegg",
               col.names = T,
               row.names = T,
               append = T)
  },
  error = function(e){
    message("nothing")
  }
)

## 2. GO
tmplist <- lapply(c("BP","CC","MF"), function(ont){
  cat(ont, "\n")
  x = enrichGO(gene = entrezid,
               OrgDb = ifelse(species=="mouse", "org.Mm.eg.db", "org.Hs.eg.db"),
               ont = ont,
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.1,
               minGSSize = 10,
               maxGSSize = 300,
               readable = TRUE)
  y = x %>% data.frame
  y$ONTOLOGY = ont
  # y = x[which(x$Count>10),] %>% data.frame
  gr = stringr::str_split_fixed(y$GeneRatio,"/",2)
  br = stringr::str_split_fixed(y$BgRatio,"/",2)
  fe = (as.numeric(gr[,1])/as.numeric(gr[,2])) / (as.numeric(br[,1])/as.numeric(br[,2]))
  y$fold_enrichment = fe
  y$generatio = as.numeric(gr[,1])/as.numeric(gr[,2])
  y$bgratio = as.numeric(br[,1])/as.numeric(br[,2])
  y
})
go_result = do.call(rbind, tmplist) %>% as.data.frame
save(go_result, file = paste0("enrich_go_list.rda"))

tryCatch({
    write.xlsx(go_result,
               file = paste0("enrich_go_combine.xlsx"),
               sheetName = "go",
               col.names = T,
               row.names = T,
               append = T)
  },
  error = function(e){
    message("nothing")
  }
)

## 4. reactome
tryCatch(
  {
    x = enrichPathway( gene = entrezid, 
                       organism = species,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.1,
                       minGSSize = 10,
                       maxGSSize = 300,
                       readable=TRUE)
    y = x %>% data.frame
    # y = x[which(x$Count>10),] %>% data.frame
    gr = stringr::str_split_fixed(y$GeneRatio,"/",2)
    br = stringr::str_split_fixed(y$BgRatio,"/",2)
    fe = (as.numeric(gr[,1])/as.numeric(gr[,2])) / (as.numeric(br[,1])/as.numeric(br[,2]))
    y$fold_enrichment = fe
    y$generatio = as.numeric(gr[,1])/as.numeric(gr[,2])
    y$bgratio = as.numeric(br[,1])/as.numeric(br[,2])
    y
  },
  error = function(e){
    NULL
  }
)
reactome_result = y
save(reactome_result, file = paste0("enrich_reactome_list.rda"))

tryCatch(
  {
    write.xlsx(reactome_result,
               file = paste0("enrich_reactome_combine.xlsx"),
               sheetName = "reactome",
               col.names = T,
               row.names = T,
               append = T)
  },
  error = function(e){
    message("nothing")
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                              
#                              enrichment plot
#                                                                              
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

enrich_label = "reactome"
q.value.filter = 0.05
df = get(paste0(enrich_label,"_result"))
df = df[which(df$qvalue < q.value.filter),]
df = df[order(df$fold_enrichment),]
df = df[!duplicated(df$Description),]
df$Description = factor(df$Description, levels = df$Description )
p <- ggplot(df, aes(x = Description, y = fold_enrichment)) +
  geom_col(aes(fill = qvalue), width = 0.5) +
  scale_fill_gradient(high = 'darkred', low = 'darkgreen') +
  # scale_fill_gradient2(colors = c('#55B047', '#FBA304', '#FF1900')) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  # theme(legend.position = 'none') + 
  coord_flip() +
  labs(x = '', y = 'Enrichment Score')
if(enrich_label=="go"){
  p = p + facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y') 
}

mis = length(which(df$qvalue < q.value.filter))
pdf(paste0("enrich_",enrich_label,"_qfilter",q.value.filter,".pdf"), 
    width = 10, height = mis * 0.15)
p
dev.off()
