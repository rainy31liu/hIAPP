library(GeneOverlap)
library(Biobase)
library(annotate)
library(biomaRt)
library(tibble)
library(rlist)
library(stringr)
library(corrplot)

setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab/T2D/")

#########get up- and down- regulated gene list
hIAPP.6wk <- read.csv("./DEG_wk6/hIAPP_results.csv", header = TRUE)
hIAPP.9wk <- read.csv("./DEG_wk6/hIAPP_results.csv", header = TRUE)

up.de.6wk <- hIAPP.6wk[hIAPP.6wk$log2FoldChange>0 & hIAPP.6wk$padj < 0.05, c(1,3,7)]
down.de.6wk <- hIAPP.6wk[hIAPP.6wk$log2FoldChange<0 & hIAPP.6wk$padj < 0.05, c(1,3,7)]
up.de.9wk <- hIAPP.9wk[hIAPP.9wk$log2FoldChange>0 & hIAPP.9wk$padj < 0.05, c(1,3,7)]
down.de.9wk <- hIAPP.9wk[hIAPP.9wk$log2FoldChange<0 & hIAPP.9wk$padj < 0.05, c(1,3,7)]

deg.lists <- tibble::tibble(degset = list(c('up.de.6wk','down.de.6wk','up.de.9wk','down.de.9wk')), 
                            genes = list(up.de.6wk$X, down.de.9wk$X, up.de.9wk$X, down.de.9wk$X))
str(deg.lists)

deg.lists <- tibble::tibble(genes = list(up.de.6wk$X, down.de.9wk$X, up.de.9wk$X, down.de.9wk$X))


###########read in markers and filter each type to individuals 
setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab/T2D/Cell_Enrichment")
cellmarkers<-read.delim("cellmarkers.txt")

peri.m <- dplyr::filter(cellmarkers, `cell.type`=="Pericytes")
mast.m <- dplyr::filter(cellmarkers, `cell.type`=="Mast cells")
psc.m <- dplyr::filter(cellmarkers, `cell.type`=="Pancreatic stellate cells")
macro.m <- dplyr::filter(cellmarkers, `cell.type`=="Macrophages")
beta.m <- dplyr::filter(cellmarkers, `cell.type`=="Beta cells")
sst.m <- dplyr::filter(cellmarkers, `cell.type`=="Delta cells")
alpha.m <- dplyr::filter(cellmarkers, `cell.type`=="Alpha cells")
ppro.m <- dplyr::filter(cellmarkers, `cell.type`=="Pancreatic progenitor cells")
pis.m <- dplyr::filter(cellmarkers, `cell.type`=="Peri-islet Schwann cells")
duct.m <- dplyr::filter(cellmarkers, `cell.type`=="Ductal cells")
# aci.m <- dplyr::filter(cellmarkers, `cell.type`=="Acinar cells")
epsi.m <- dplyr::filter(cellmarkers, `cell.type`=="Epsilon cells")
pp.m <- dplyr::filter(cellmarkers, `cell.type`=="Gamma (PP) cells")
vec.m <- dplyr::filter(cellmarkers, `cell.type`=="Endothelial cells")

alpha.list <- alpha.m$official.gene.symbol
beta.list <- beta.m$official.gene.symbol
sst.list <- sst.m$official.gene.symbol
pp.list <- pp.m$official.gene.symbol
epsi.list <- epsi.m$official.gene.symbol
ppro.list <- ppro.m$official.gene.symbol
duct.list <- duct.m$official.gene.symbol
# aci.list <- aci.m$official.gene.symbol
psc.list <- psc.m$official.gene.symbol
pis.list <- pis.m$official.gene.symbol
macro.list <- macro.m$official.gene.symbol
peri.list <- peri.m$official.gene.symbol
vec.list <- vec.m$official.gene.symbol


######### create list
cell.marker.list <- list(alpha.list, beta.list, sst.list, pp.list, epsi.list, ppro.list, 
                         duct.list, psc.list, pis.list, macro.list, peri.list, vec.list)
names(cell.marker.list) <- c("Alpha","Beta","Delta","Gamma (PP)","Epsilon","Pancreatic progenitor",
                             "Ductal","Stellate","Schwann","Macrophage","Pericyte","Endothelial")

########## gene list for each contrast
up.de.6wk.list <- toupper(up.de.6wk$X)
down.de.6wk.list <- toupper(down.de.6wk$X)
up.de.9wk.list <- toupper(up.de.9wk$X)
down.de.9wk.list <- toupper(down.de.9wk$X)


deg.list <- list(up.de.6wk.list,down.de.6wk.list,up.de.9wk.list,down.de.9wk.list)
names(deg.list) <- c("Up 6wk hIAPP vs. WT","Down 6wk hIAPP vs. WT","Up 9wk hIAPP vs. WT",
                     "Down 9wk hIAPP vs. WT")

deg.cell.GOMobj <- newGOM(gsetA = deg.list, gsetB = cell.marker.list)


########### get odds ratio and p-values matrix 
deg.cell.GOMobj <- newGOM(gsetA = deg.list, gsetB = cell.marker.list, spec = "mm9.gene")
deg.cell.ormat <- getMatrix(deg.cell.GOMobj, "odds.ratio")
deg.cell.pmat <- getMatrix(deg.cell.GOMobj, "pval")
or.as.pmat <- as_tibble(deg.cell.ormat)

# deg.cell.pmat.adj <- deg.cell.pmat %>% as.vector() %>% p.adjust(method = "fdr") %>% matrix(ncol = 13)
deg.cell.pmat.adj <- deg.cell.pmat %>% as.vector() %>% p.adjust(method = "fdr") %>% matrix(ncol = 12)
rownames(deg.cell.pmat.adj) <- rownames(deg.cell.pmat)
colnames(deg.cell.pmat.adj) <- colnames(deg.cell.pmat)
or.as.pmat <- as_tibble(deg.cell.ormat)

or.as.pmat$Alpha <- c(NA, 3.99,4.23,13.81)
or.as.pmat$Beta <- c(NA,7.05,NA,32.35)
or.as.pmat$Delta <- c(NA,6.33,NA,9.28)
or.as.pmat$`Gamma (PP)` <- c(NA,4.06,NA,NA)
or.as.pmat$Epsilon <-c(4.13,NA,NA,NA)
or.as.pmat$`Pancreatic progenitor` <- NA
or.as.pmat$Ductal <- c(NA,4.16,NA,NA)
or.as.pmat$Stellate <- c(6.97,NA,26.05,NA)
or.as.pmat$Schwann <- c(NA,NA,8.53,NA)
or.as.pmat$Macrophage <- c(8.86,NA,12.66,NA)
or.as.pmat$Pericyte <- c(9.39,NA,17.48,NA)
or.as.pmat$Endothelial <- c(NA,1.99,4.53,NA)
rownames(or.as.pmat) <- c("Up 6wk hIAPP vs. WT","Down 6wk hIAPP vs. WT","Up 9wk hIAPP vs. WT",
                          "Down 9wk hIAPP vs. WT")

tiff("6wk_hIAPP_noAcinar_1.tiff", width = 12, height = 5, units = "in", res = 720)
corrplot(-log10(deg.cell.pmat.adj),  order = c("original"), cl.offset = 6, mar = c(1,2,2,2), 
         addgrid.col = "black", method = "color", is.corr = F, col=colorRampPalette(c("white","firebrick3"))(200), 
         cl.lim = NULL, cl.length = 3, tl.col = "black", tl.srt = 45, tl.cex = 1, cl.cex = 1, 
         p.mat =  matrix(or.as.pmat),
         sig.level = 1, insig = "p-value", number.cex = 0.5, number.digits = 2)
dev.off()

