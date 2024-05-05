library(RRHO2)
library(lattice)
library(gridExtra)

setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab/T2D")
hIAPP.6wk <- read.csv("./DEG_wk6/hIAPP_results.csv")
hIAPP.9wk <- read.csv("./DEG_wk9/hIAPP_results.csv")
T2D.diab.ctrl <- read.csv("./DEG_T2D/diab_results.csv")
T2D.pre.ctrl <- read.csv("./DEG_T2D/pre_results.csv")

hIAPP.6wk$metric <- -log10(hIAPP.6wk$padj) * sign(hIAPP.6wk$log2FoldChange)
hIAPP.9wk$metric <- -log10(hIAPP.9wk$padj) * sign(hIAPP.9wk$log2FoldChange)
T2D.diab.ctrl$metric <- -log10(T2D.diab.ctrl$padj) * sign(T2D.diab.ctrl$log2FoldChange)
T2D.pre.ctrl$metric <- -log10(T2D.pre.ctrl$padj) * sign(T2D.pre.ctrl$log2FoldChange)

hIAPP.6wk <- hIAPP.6wk[!is.na(hIAPP.6wk$padj),]
hIAPP.9wk <- hIAPP.9wk[!is.na(hIAPP.9wk$padj),]
T2D.diab.ctrl <- T2D.diab.ctrl[!is.na(T2D.diab.ctrl$padj),]
T2D.pre.ctrl <- T2D.pre.ctrl[!is.na(T2D.pre.ctrl$padj),]


jet.colors  <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

###################### hIAPP vs WT 6wk, hIAPP vs WT 9wk#########################
hIAPP.filtered <- merge(hIAPP.6wk, hIAPP.9wk, by.x = "X", by.y = "X")
hIAPP.filtered <- hIAPP.filtered[!duplicated(hIAPP.filtered$X),]

rrho.hIAPP <- RRHO2_initialize(list1 = as.data.frame(hIAPP.filtered[,c(1,8)]),
                               list2 = as.data.frame(hIAPP.filtered[,c(1,15)]),
                               labels = c("hIAPP vs. WT 6wk", "hIAPP vs. WT 9wk"),log10.ind=TRUE)
rrho.hIAPP$hypermat[rrho.hIAPP$hypermat>100] <- 100

tiff("./RRHO/hIAPP_9wk_6wk.tiff", bg = "transparent", res = 600, width = 3.5, height = 3.5, units = "in", compression = "lzw")
levelplot(rrho.hIAPP$hypermat, xlab="hIAPP vs. WT 6wk", ylab="hIAPP vs. WT 9wk", col.regions=jet.colors, 
          at=seq(0,100), pretty=TRUE, scales = list(tck = c(1,0)))
dev.off()

## conversion ##################################################################

load("../gene_symbol_conversions.rda")
module_list <- list("wk6" = hIAPP.6wk, "wk9" = hIAPP.9wk)

for (df_name in names(module_list)) {
  df <- module_list[[df_name]]
  colnames <- colnames(df)
  df <- left_join(df, Mouse_Human, by = c("X" = "mouse_symbol"))
  df <- df[complete.cases(df), ]
  df <- select(df, -X, human_symbol)
  df <- df[, c(8,1:7)]
  colnames(df) <- colnames
  module_list[[df_name]] <- df
}
hIAPP.6wk <- module_list[[1]]
hIAPP.9wk <- module_list[[2]]

###################### hIAPP vs WT 6wk, Prediabetics vs. Control################
hIAPP.pre.filtered <- merge(hIAPP.6wk, T2D.pre.ctrl, by.x = "X", by.y = "X")
hIAPP.pre.filtered <- hIAPP.pre.filtered[!duplicated(hIAPP.pre.filtered$X),]

rrho.hIAPP.pre <- RRHO2_initialize(list1 = as.data.frame(hIAPP.pre.filtered[,c(1,8)]),
                               list2 = as.data.frame(hIAPP.pre.filtered[,c(1,15)]),
                               labels = c("hIAPP vs. WT 6wk", "Prediabetics vs. Control"),log10.ind=TRUE)
rrho.hIAPP.pre$hypermat[rrho.hIAPP.pre$hypermat>70] <- 70

tiff("./RRHO/hIAPP_6wk_Prediabetics.tiff", bg = "transparent", res = 600, width = 3.5, height = 3.5, units = "in", compression = "lzw")
levelplot(rrho.hIAPP.pre$hypermat, xlab="hIAPP vs. WT 6wk", ylab="Prediabetics vs. Control", col.regions=jet.colors, at=seq(0,70), pretty=TRUE, scales = list(tck = c(1,0)))
dev.off()

###################### hIAPP vs WT 6wk, Diabetics vs. Control###################
hIAPP.diab.filtered <- merge(hIAPP.6wk, T2D.diab.ctrl, by.x = "X", by.y = "X")
hIAPP.diab.filtered <- hIAPP.diab.filtered[!duplicated(hIAPP.diab.filtered$X),]

rrho.hIAPP.diab <- RRHO2_initialize(list1 = as.data.frame(hIAPP.diab.filtered[,c(1,8)]),
                                   list2 = as.data.frame(hIAPP.diab.filtered[,c(1,15)]),
                                   labels = c("hIAPP vs. WT 6wk", "T2D vs. Control"),log10.ind=TRUE)
rrho.hIAPP.diab$hypermat[rrho.hIAPP.diab$hypermat>70] <- 70

tiff("./RRHO/hIAPP_6wk_T2D.tiff", bg = "transparent", res = 600, width = 3.5, height = 3.5, units = "in", compression = "lzw")
levelplot(rrho.hIAPP.diab$hypermat, xlab="hIAPP vs. WT 6wk", ylab="T2D vs. Control", col.regions=jet.colors, at=seq(0,70), pretty=TRUE, scales = list(tck = c(1,0)))
dev.off()
