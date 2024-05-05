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
raw.counts.wk6 <- read_csv("txi_hIAPP_wk6.csv") %>% 
  data.frame %>% filter(!grepl("Mar|Sep|Rik|LOC", gene))

########### keep counts columns only
rownames(raw.counts.wk6) <- raw.counts.wk6$gene
counts.only.wk6 <- dplyr::select(raw.counts.wk6, -gene) %>% as.matrix 

########### filter counts and standard deviation estimates
counts.only.wk6 <- counts.only.wk6[rowSums(counts.only.wk6) > 1, ]
counts.mad.wk6 <- rowMads(counts.only.wk6)
counts.only.wk6 <- counts.only.wk6[counts.mad.wk6 > 0,]

########## create DGEList object
counts.dge.wk6 <- DGEList(counts.only.wk6)
counts.dgenorm.wk6 <- calcNormFactors(counts.dge.wk6)

########## use column names to create pheno data
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
contrasts.matrix.wk6 <- model.matrix( ~ 0 + Combined, colnames.split.wk6)
colnames(contrasts.matrix.wk6) %<>% str_replace('Combined', "")
contrasts.all.wk6 <- makeContrasts("hIAPP_vs_WT" = hIAPP_WT - WT_WT, 
                                   "rIAPP_vs_WT" = rIAPP_WT - WT_WT, 
                                   "WT_vs_WT" = WT_WT - WT_WT,
                                   levels = contrasts.matrix.wk6)

########### use limma voom to fit model and save object for later use
voom.tmm.wk6 <- voom(counts.dgenorm.wk6, design = contrasts.matrix.wk6, normalize = "quantile")

##############ratio
voom.tmm.expr.wk6<-voom.tmm.wk6$E
all.samples.wk6<-as.data.frame(voom.tmm.expr.wk6)

######get average for each group
wt.wk6<-all.samples.wk6[,colnames.split.wk6$Combined == "WT_WT"]
wtM.wk6 <-rowMeans(wt.wk6)
rIAPP.wk6 <- all.samples.wk6[,colnames.split.wk6$Combined== "rIAPP_WT" ] 
hIAPP.wk6<-all.samples.wk6[,colnames.split.wk6$Combined== "hIAPP_WT" ]

hIAPP.wt.wk6<- hIAPP.wk6- wtM.wk6
colnames(hIAPP.wt.wk6)<-paste(colnames(hIAPP.wt.wk6), "_vs_WT" ,sep="")
rIAPP.wt.wk6<- rIAPP.wk6- wtM.wk6
colnames(rIAPP.wt.wk6)<-paste(colnames(rIAPP.wt.wk6), "_vs_WT" ,sep="")

ratio.expr.wk6<-cbind(Symbol=rownames(voom.tmm.expr.wk6), hIAPP.wt.wk6, rIAPP.wt.wk6)
# write.csv(ratio.expr.wk6, "ratio_wk6.csv")

###### generate graph

### cholosterol genes
lipid_df <- read.csv("./pathway_genes/lipid.csv")
lipid <- unique(lipid_df$hIAPP)
lipid <- lipid[!is.na(lipid) & lipid != ""]
lipid<-as.data.frame(lipid)
colnames(lipid)<-"Symbol"

###### ratios from mouse
lipid.expr.wk6 <- filter(ratio.expr.wk6, Symbol %in% lipid$Symbol)
lipid.expr.wk6 <- lipid.expr.wk6[match(lipid$Symbol, lipid.expr.wk6$Symbol), ]
lipid.expr <- cbind(lipid.expr.wk6[,-1])


labeltiapp2dimp<-c(rep("",2),"       hIAPP_6wk",rep("",6),"rIAPP_6wk       ",rep("",2))
heatmap_width <- 240
heatmap_height <- 400
cell_width <- 20
cell_height <- 20


###### create graph
### lipid
heatmap_width <- 240
heatmap_height <- 420
cell_width <- 20
cell_height <- 20
cell_height <- heatmap_height/nrow(lipid.expr)
lipid.myrange <- max(abs(lipid.expr));
lipid.breaks = seq(-2.5, 2.5, length.out = 30)
lipid.heatmap <- lipid.expr[order(rowMeans(lipid.expr[, 1:6])), ]
pdf("./wk6_pathway/lipid_6wk.pdf", width = 10, height = 10)
pheatmap(lipid.heatmap, main = "Cholesterol",
         breaks = lipid.breaks,legend = T,  gaps_col = 6,labels_col = labeltiapp2dimp, 
         labels_row=toupper(rownames(lipid.heatmap)),cluster_cols = F, color = colorRampPalette(c("blue","white","red"))(30), 
         fontsize_col = 10, angle_col = 0, border_color = NA, scale = "none", cluster_rows = F,keysize = 0.2,
         cellwidth = cell_width, cellheight = cell_height, width = heatmap_width, height = heatmap_height)
dev.off()


