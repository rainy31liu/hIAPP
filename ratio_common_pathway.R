library(magrittr)
library(readr)
library(dplyr)
library(limma)
library(edgeR)
library(matrixStats)
library(stringr)
library(pheatmap)

########### read in count matrix
setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab/T2D/ratio")
raw.counts.wk9 <- read_csv("txi_hIAPP_wk9.csv") %>% 
  data.frame %>% filter(!grepl("Mar|Sep|Rik|LOC", gene))
raw.counts.wk6 <- read_csv("txi_hIAPP_wk6.csv") %>% 
  data.frame %>% filter(!grepl("Mar|Sep|Rik|LOC", gene))

########### keep counts columns only
rownames(raw.counts.wk9) <- raw.counts.wk9$gene
counts.only.wk9 <- dplyr::select(raw.counts.wk9, -gene) %>% as.matrix 

rownames(raw.counts.wk6) <- raw.counts.wk6$gene
counts.only.wk6 <- dplyr::select(raw.counts.wk6, -gene) %>% as.matrix 

########### filter counts and standard deviation estimates
counts.only.wk9 <- counts.only.wk9[rowSums(counts.only.wk9) > 1, ]
counts.mad.wk9 <- rowMads(counts.only.wk9)
counts.only.wk9 <- counts.only.wk9[counts.mad.wk9 > 0,]

counts.only.wk6 <- counts.only.wk6[rowSums(counts.only.wk6) > 1, ]
counts.mad.wk6 <- rowMads(counts.only.wk6)
counts.only.wk6 <- counts.only.wk6[counts.mad.wk6 > 0,]

########## create DGEList object
counts.dge.wk9 <- DGEList(counts.only.wk9)
counts.dgenorm.wk9 <- calcNormFactors(counts.dge.wk9)

counts.dge.wk6 <- DGEList(counts.only.wk6)
counts.dgenorm.wk6 <- calcNormFactors(counts.dge.wk6)

########## use column names to create pheno data
colnames.split.wk9 <- str_split_fixed(colnames(counts.only.wk9), "_", 3) %>% as_tibble
colnames.split.wk9

colnames.split.wk9$Fix <- colnames.split.wk9$V3 == ""
colnames.split.wk9[colnames.split.wk9$Fix,]$V3 <- colnames.split.wk9[colnames.split.wk9$Fix,]$V2
colnames.split.wk9[colnames.split.wk9$Fix,]$V2 <- "WT"
colnames.new.wk9 <- str_c(colnames.split.wk9$V1, colnames.split.wk9$V2, colnames.split.wk9$V3, sep = "_")
colnames(counts.dgenorm.wk9) <- colnames.new.wk9
colnames(colnames.split.wk9) <- c("IAPP", "WT", "Id", "Fix")
colnames.split.wk9$IAPP %<>% factor(levels = c("WT", "hIAPP", "rIAPP"))
colnames.split.wk9$WT %<>% factor(levels = c("WT"))
colnames.split.wk9$Combined <- str_c(colnames.split.wk9$IAPP, colnames.split.wk9$WT, sep = "_")

colnames.split.wk6 <- str_split_fixed(colnames(counts.only.wk6), "_", 3) %>% as_tibble
colnames.split.wk6

colnames.split.wk6$Fix <- colnames.split.wk6$V3 == ""
colnames.split.wk6[colnames.split.wk6$Fix,]$V3 <- colnames.split.wk6[colnames.split.wk6$Fix,]$V2
colnames.split.wk6[colnames.split.wk6$Fix,]$V2 <- "WT"
colnames.new.wk6 <- str_c(colnames.split.wk6$V1, colnames.split.wk6$V2, colnames.split.wk6$V3, sep = "_")
colnames(counts.dgenorm.wk6) <- colnames.new.wk6
colnames(colnames.split.wk6) <- c("IAPP", "WT", "Id", "Fix")
colnames.split.wk6$IAPP %<>% factor(levels = c("WT", "hIAPP", "rIAPP"))
colnames.split.wk6$WT %<>% factor(levels = c("WT"))
colnames.split.wk6$Combined <- str_c(colnames.split.wk6$IAPP, colnames.split.wk6$WT, sep = "_")

########### DEG analysis
contrasts.matrix.wk9 <- model.matrix( ~ 0 + Combined, colnames.split.wk9)
colnames(contrasts.matrix.wk9) %<>% str_replace('Combined', "")
contrasts.all.wk9 <- makeContrasts("hIAPP_vs_WT" = hIAPP_WT - WT_WT, 
                                   "rIAPP_vs_WT" = rIAPP_WT - WT_WT, 
                                   "WT_vs_WT" = WT_WT - WT_WT,
                                   levels = contrasts.matrix.wk9)

contrasts.matrix.wk6 <- model.matrix( ~ 0 + Combined, colnames.split.wk6)
colnames(contrasts.matrix.wk6) %<>% str_replace('Combined', "")
contrasts.all.wk6 <- makeContrasts("hIAPP_vs_WT" = hIAPP_WT - WT_WT, 
                                   "rIAPP_vs_WT" = rIAPP_WT - WT_WT, 
                                   "WT_vs_WT" = WT_WT - WT_WT,
                                   levels = contrasts.matrix.wk6)

########### use limma voom to fit model and save object for later use
voom.tmm.wk9 <- voom(counts.dgenorm.wk9, design = contrasts.matrix.wk9, normalize = "quantile")
voom.tmm.wk6 <- voom(counts.dgenorm.wk6, design = contrasts.matrix.wk6, normalize = "quantile")

##############ratio
voom.tmm.expr.wk9<-voom.tmm.wk9$E
all.samples.wk9<-as.data.frame(voom.tmm.expr.wk9)

voom.tmm.expr.wk6<-voom.tmm.wk6$E
all.samples.wk6<-as.data.frame(voom.tmm.expr.wk6)

######get average for each group
wt.wk9<-all.samples.wk9[,colnames.split.wk9$Combined == "WT_WT"]
wtM.wk9 <-rowMeans(wt.wk9)
rIAPP.wk9 <- all.samples.wk9[,colnames.split.wk9$Combined== "rIAPP_WT" ] 
hIAPP.wk9<-all.samples.wk9[,colnames.split.wk9$Combined== "hIAPP_WT" ]

hIAPP.wt.wk9<- hIAPP.wk9- wtM.wk9
colnames(hIAPP.wt.wk9)<-paste(colnames(hIAPP.wt.wk9), "_vs_WT" ,sep="")
rIAPP.wt.wk9<- rIAPP.wk9- wtM.wk9
colnames(rIAPP.wt.wk9)<-paste(colnames(rIAPP.wt.wk9), "_vs_WT" ,sep="")

ratio.expr.wk9<-cbind(Symbol=rownames(voom.tmm.expr.wk9), hIAPP.wt.wk9, rIAPP.wt.wk9)
write.csv(ratio.expr.wk9, "ratio_hIAPP_wk9.csv")

wt.wk6<-all.samples.wk6[,colnames.split.wk6$Combined == "WT_WT"]
wtM.wk6 <-rowMeans(wt.wk6)
rIAPP.wk6 <- all.samples.wk6[,colnames.split.wk6$Combined== "rIAPP_WT" ] 
hIAPP.wk6<-all.samples.wk6[,colnames.split.wk6$Combined== "hIAPP_WT" ]

hIAPP.wt.wk6<- hIAPP.wk6- wtM.wk6
colnames(hIAPP.wt.wk6)<-paste(colnames(hIAPP.wt.wk6), "_vs_WT" ,sep="")
rIAPP.wt.wk6<- rIAPP.wk6- wtM.wk6
colnames(rIAPP.wt.wk6)<-paste(colnames(rIAPP.wt.wk6), "_vs_WT" ,sep="")

ratio.expr.wk6<-cbind(Symbol=rownames(voom.tmm.expr.wk6), hIAPP.wt.wk6, rIAPP.wt.wk6)
write.csv(ratio.expr.wk6, "ratio_hIAPP_wk6.csv")

###### generate graph

### Unfolded protein response genes
upr<-c("Ppp1r15a","Stc2","Xbp1","Pdia6","Bhlha15","Cdk5rap3","Creb3","Dnajc3","Edem3","Eif2ak3","Manf","Syvn1","Dab2ip")
### Cell cycle genes
ccy<-c("Nek2","Cenpe","Anln","Hmmr","Kif2c","Aurkb","Cep55","Bub1","Cdk1","Ccnb1","Bub1b","Spag5","G2e3","Plk1","Ccnb2","Aspm","Mcm2")
### inflammation genes
inf<-c("Cx3cr1","Lgals3","Axl","Cd84","Slc11a1","Cd4","Cd68","Clec7a","Trem2","Csf1r","Cxcl10","Fcgr2b","Itgam","Itgax","Mpeg1","Tlr7","Tnfrsf1b","Irf8","Ncf1","Socs3","Tnfrsf1a","Cxcl12")
### differentiation genes
diffe<-c("Ezh1","Jph3","Igf2","Neurod1","Nkx6-1","Pcsk2","Pdx1","Pfkfb2","Rgs16","Slc2a2","Slc30a8","Npy","G6pc2","Pcsk1","Edaradd","Gjd2","Casr","Nkx2-2","Efna5","Mafa","Hopx")
### cholosterol genes
lipid<-c("Acly","Dhcr24","Fdft1","Fgf1","Hmgcr","Hmgcs1","Hsd17b7","Insig1","Lbr","Lss","Msmo1","Mvd","Mvk","Nr0b2","Nsdhl","Osbpl5","Prkaa1","Sqle","Srebf1","Srebf2","Stard4")

upr<-as.data.frame(upr)
colnames(upr)<-"Symbol"
ccy<-as.data.frame(ccy)
colnames(ccy)<-"Symbol"
inf<-as.data.frame(inf)
colnames(inf)<-"Symbol"
diffe<-as.data.frame(diffe)
colnames(diffe)<-"Symbol"
lipid<-as.data.frame(lipid)
colnames(lipid)<-"Symbol"

###### ratios from mouse
dedif.expr.wk9 <- filter(ratio.expr.wk9, Symbol %in% diffe$Symbol)
dedif.expr.wk9 <- dedif.expr.wk9[match(diffe$Symbol, dedif.expr.wk9$Symbol), ]
dedif.expr.wk6 <- filter(ratio.expr.wk6, Symbol %in% diffe$Symbol)
dedif.expr.wk6 <- dedif.expr.wk6[match(diffe$Symbol, dedif.expr.wk6$Symbol), ]
dedif.expr <- cbind(dedif.expr.wk6[,-1],dedif.expr.wk9[,-1])

cellcyc.expr.wk9 <- filter(ratio.expr.wk9, Symbol %in% ccy$Symbol)
cellcyc.expr.wk9 <- cellcyc.expr.wk9[match(ccy$Symbol, cellcyc.expr.wk9$Symbol), ]
cellcyc.expr.wk6 <- filter(ratio.expr.wk6, Symbol %in% ccy$Symbol)
cellcyc.expr.wk6 <- cellcyc.expr.wk6[match(ccy$Symbol, cellcyc.expr.wk6$Symbol), ]
cellcyc.expr <- cbind(cellcyc.expr.wk6[,-1],cellcyc.expr.wk9[,-1])

inf.expr.wk9 <- filter(ratio.expr.wk9, Symbol %in% inf$Symbol)
inf.expr.wk9 <- inf.expr.wk9[match(inf$Symbol, inf.expr.wk9$Symbol), ]
inf.expr.wk6 <- filter(ratio.expr.wk6, Symbol %in% inf$Symbol)
inf.expr.wk6 <- inf.expr.wk6[match(inf$Symbol, inf.expr.wk6$Symbol), ]
inf.expr <- cbind(inf.expr.wk6[,-1],inf.expr.wk9[,-1])

upr.expr.wk9 <- filter(ratio.expr.wk9, Symbol %in% upr$Symbol)
upr.expr.wk9 <- upr.expr.wk9[match(upr$Symbol, upr.expr.wk9$Symbol), ]
upr.expr.wk6 <- filter(ratio.expr.wk6, Symbol %in% upr$Symbol)
upr.expr.wk6 <- upr.expr.wk6[match(upr$Symbol, upr.expr.wk6$Symbol), ]
upr.expr <- cbind(upr.expr.wk6[,-1],upr.expr.wk9[,-1])

lipid.expr.wk9 <- filter(ratio.expr.wk9, Symbol %in% lipid$Symbol)
lipid.expr.wk9 <- lipid.expr.wk9[match(lipid$Symbol, lipid.expr.wk9$Symbol), ]
lipid.expr.wk6 <- filter(ratio.expr.wk6, Symbol %in% lipid$Symbol)
lipid.expr.wk6 <- lipid.expr.wk6[match(lipid$Symbol, lipid.expr.wk6$Symbol), ]
lipid.expr <- cbind(lipid.expr.wk6[,-1],lipid.expr.wk9[,-1])

labeltiapp2dimp<-c(rep("",2),"       hIAPP_6wk",rep("",3),"hIAPP_9wk","")
heatmap_width <- 180
heatmap_height <- 260
cell_width <- 20

###### create graph

### upr
cell_height <- heatmap_height/nrow(upr.expr)
upr.myrange <- max(abs(upr.expr));
upr.breaks = seq(-1.5, 1.5, length.out = 30)
upr.heatmap <- upr.expr[,c(1:6,13:15)]
upr.heatmap <- upr.heatmap[order(-rowMeans(upr.heatmap[, 1:6])), ]
pdf("./hIAPP_pathway/UPR_hIAPP.pdf", width = 5, height = 5)
pheatmap(upr.heatmap, main = "Unfolded protein response",
         breaks = upr.breaks,legend = T,  gaps_col = 6,labels_col = labeltiapp2dimp, 
         labels_row=toupper(rownames(upr.heatmap)),cluster_cols = F, color = colorRampPalette(c("blue","white","red"))(30), 
         fontsize_col = 10, angle_col = 0, border_color = NA, scale = "none", cluster_rows = F,keysize = 0.2,
         cellwidth = cell_width, cellheight = cell_height, width = heatmap_width, height = heatmap_height,)
dev.off()

### cell cycle
cell_height <- heatmap_height/nrow(cellcyc.expr)
ccy.myrange <- max(abs(cellcyc.expr));
ccy.breaks = seq(-2.5, 2.5, length.out = 30)
ccy.heatmap <- cellcyc.expr[,c(1:6,13:15)]
ccy.heatmap <- ccy.heatmap[order(-rowMeans(ccy.heatmap[, c(1:4, 7:9)])), ]
pdf("./hIAPP_pathway/Cell_cycle_hIAPP.pdf", width = 5, height = 5)
pheatmap(ccy.heatmap, main = "Cell Cycle",
         breaks = ccy.breaks,legend = T,  gaps_col = 6,labels_col = labeltiapp2dimp, 
         labels_row=toupper(rownames(ccy.heatmap)),cluster_cols = F, color = colorRampPalette(c("blue","white","red"))(30), 
         fontsize_col = 10, angle_col = 0, border_color = NA, scale = "none", cluster_rows = F,keysize = 0.2,
         cellwidth = cell_width, cellheight = cell_height, width = heatmap_width, height = heatmap_height,)
dev.off()

### differentiation
cell_height <- heatmap_height/nrow(dedif.expr)
ddf.myrange <- max(abs(dedif.expr));
ddf.breaks = seq(-3.2, 3.2, length.out = 30)
ddf.heatmap <- dedif.expr[,c(1:6,13:15)]
ddf.heatmap <- ddf.heatmap[order(rowMeans(ddf.heatmap[, 1:6])), ]
pdf("./hIAPP_pathway/Dedifferentiation_hIAPP.pdf", width = 5, height = 5)
pheatmap(ddf.heatmap, main = "Differentiation",
         breaks = ddf.breaks,legend = T,  gaps_col = 6,labels_col = labeltiapp2dimp, 
         labels_row=toupper(rownames(ddf.heatmap)),cluster_cols = F, color = colorRampPalette(c("blue","white","red"))(30), 
         fontsize_col = 10, angle_col = 0, border_color = NA, scale = "none", cluster_rows = F,keysize = 0.2,
         cellwidth = cell_width, cellheight = cell_height, width = heatmap_width, height = heatmap_height,)
dev.off()

### inflammation
cell_height <- heatmap_height/nrow(inf.expr)
inf.myrange <- max(abs(inf.expr));
inf.breaks = seq(-4.5, 4.5, length.out = 30)
inf.heatmap <- inf.expr[,c(1:6,13:15)]
inf.heatmap <- inf.heatmap[order(-rowMeans(inf.heatmap[, 1:6])), ]
pdf("./hIAPP_pathway/Inflammation_hIAPP.pdf", width = 5, height = 5)
pheatmap(inf.heatmap, main = "Inflammation",
         breaks = inf.breaks,legend = T,  gaps_col = 6,labels_col = labeltiapp2dimp, 
         labels_row=toupper(rownames(inf.heatmap)),cluster_cols = F, color = colorRampPalette(c("blue","white","red"))(30), 
         fontsize_col = 10, angle_col = 0, border_color = NA, scale = "none", cluster_rows = F,keysize = 0.2,
         cellwidth = cell_width, cellheight = cell_height, width = heatmap_width, height = heatmap_height,)
dev.off()

### lipid
cell_height <- heatmap_height/nrow(lipid.expr)
lipid.myrange <- max(abs(lipid.expr))
lipid.breaks = seq(-1.5, 1.5, length.out = 30)
pdf("./hIAPP_pathway/cholosterol_hIAPP.pdf", width = 5, height = 5)
pheatmap(lipid.expr[,c(1:6,13:15)], main = "Cholesterol",
         breaks = lipid.breaks,legend = T,  gaps_col = 6,labels_col = labeltiapp2dimp, 
         labels_row=toupper(rownames(lipid.expr)),cluster_cols = F, color = colorRampPalette(c("blue","white","red"))(30), 
         fontsize_col = 10, angle_col = 0, border_color = NA, scale = "none", cluster_rows = F,keysize = 0.2,
         cellwidth = cell_width, cellheight = cell_height, width = heatmap_width, height = heatmap_height,)
dev.off()
