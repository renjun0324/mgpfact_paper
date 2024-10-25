
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

## add new label
x = ct$Cluster
x[which(x=="P60")] = "HM"
x[which(x=="P7_C0")] = "IM"
x[which(x=="P7_C1")] = "PAM"
x[which(x=="P7-Gpnmb+Clec7a+")] = "PAM"
ct@MetaData$celltype8 = x
d2_label = c("day", "truetime", "celltype", "p7_ident", "celltype6", "celltype7", "celltype8")
ct = GetMURPMapLabel(ct, labels = d2_label)

## 

## save 
save(ct, file = "ct.rda")
writeSettings(ct)