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
write.csv(ratio.expr.wk9, "ratio_wk9.csv")

wt.wk6<-all.samples.wk6[,colnames.split.wk6$Combined == "WT_WT"]
wtM.wk6 <-rowMeans(wt.wk6)
rIAPP.wk6 <- all.samples.wk6[,colnames.split.wk6$Combined== "rIAPP_WT" ] 
hIAPP.wk6<-all.samples.wk6[,colnames.split.wk6$Combined== "hIAPP_WT" ]

hIAPP.wt.wk6<- hIAPP.wk6- wtM.wk6
colnames(hIAPP.wt.wk6)<-paste(colnames(hIAPP.wt.wk6), "_vs_WT" ,sep="")
rIAPP.wt.wk6<- rIAPP.wk6- wtM.wk6
colnames(rIAPP.wt.wk6)<-paste(colnames(rIAPP.wt.wk6), "_vs_WT" ,sep="")

ratio.expr.wk6<-cbind(Symbol=rownames(voom.tmm.expr.wk6), hIAPP.wt.wk6, rIAPP.wt.wk6)
write.csv(ratio.expr.wk6, "ratio_wk6.csv")

###### generate graph
### hormone secretion
hormone_df <- read.csv("./pathway_genes/hormone_secretion.csv")
hormone <- unique(c(hormone_df$X6wk, hormone_df$X9wk))
# hormone <- unique(hormone_df$X6wk)
hormone <- hormone[!is.na(hormone) & hormone != ""]

### mRNA processing
RNA_df <- read.csv("./pathway_genes/mRNA_processing.csv")
RNA <- unique(c(RNA_df$X6wk, RNA_df$X9wk))
# RNA <- unique(RNA_df$X6wk)
RNA <- RNA[!is.na(RNA) & RNA != ""]

### ribosome biogenesis
ribo_df <- read.csv("./pathway_genes/ribosome_biogenesis.csv")
ribo <- unique(c(ribo_df$X6wk, ribo_df$X9wk))
# ribo <- unique(ribo_df$X6wk)
ribo <- ribo[!is.na(ribo) & ribo != ""]

hormone<-as.data.frame(hormone)
colnames(hormone)<-"Symbol"
RNA<-as.data.frame(RNA)
colnames(RNA)<-"Symbol"
ribo<-as.data.frame(ribo)
colnames(ribo)<-"Symbol"

###### ratios from mouse
hormone.expr.wk9 <- filter(ratio.expr.wk9, Symbol %in% hormone$Symbol)
hormone.expr.wk9 <- hormone.expr.wk9[match(hormone$Symbol, hormone.expr.wk9$Symbol), ]
hormone.expr.wk6 <- filter(ratio.expr.wk6, Symbol %in% hormone$Symbol)
hormone.expr.wk6 <- hormone.expr.wk6[match(hormone$Symbol, hormone.expr.wk6$Symbol), ]
hormone.expr <- cbind(hormone.expr.wk6[,-1],hormone.expr.wk9[,-1])

RNA.expr.wk9 <- filter(ratio.expr.wk9, Symbol %in% RNA$Symbol)
RNA.expr.wk9 <- RNA.expr.wk9[match(RNA$Symbol, RNA.expr.wk9$Symbol), ]
RNA.expr.wk6 <- filter(ratio.expr.wk6, Symbol %in% RNA$Symbol)
RNA.expr.wk6 <- RNA.expr.wk6[match(RNA$Symbol, RNA.expr.wk6$Symbol), ]
RNA.expr <- cbind(RNA.expr.wk6[,-1],RNA.expr.wk9[,-1])

ribo.expr.wk9 <- filter(ratio.expr.wk9, Symbol %in% ribo$Symbol)
ribo.expr.wk9 <- ribo.expr.wk9[match(ribo$Symbol, ribo.expr.wk9$Symbol), ]
ribo.expr.wk6 <- filter(ratio.expr.wk6, Symbol %in% ribo$Symbol)
ribo.expr.wk6 <- ribo.expr.wk6[match(ribo$Symbol, ribo.expr.wk6$Symbol), ]
ribo.expr <- cbind(ribo.expr.wk6[,-1],ribo.expr.wk9[,-1])

labeltiapp2dimp<-c(rep("",2),"       hIAPP_6wk",rep("",3),"hIAPP_9wk","")
heatmap_width <- 180
heatmap_height <- 450
cell_width <- 20

###### create graph
### hormone secretion
cell_height <- heatmap_height/nrow(hormone.expr)
hormone.myrange <- max(abs(hormone.expr));
hormone.breaks = seq(-3, 3, length.out = 30)
hormone.heatmap <- hormone.expr[,c(1:6,13:15)]
hormone.heatmap <- hormone.heatmap[order(rowMeans(hormone.heatmap[, 1:6])), ]
hormone.matrix <- matrix("", nrow = nrow(hormone.heatmap), ncol = ncol(hormone.heatmap))
pdf("./hIAPP_pathway/hormone_secretion_hIAPP.pdf", width = 8, height = 8)
pheatmap(hormone.heatmap, main = "Hormone Secretion",
         breaks = hormone.breaks, legend = T, gaps_col = 6,labels_col = labeltiapp2dimp,
         labels_row=toupper(rownames(hormone.heatmap)),cluster_cols = F, color = colorRampPalette(c("blue","white","red"))(30),
         fontsize_col = 10, fontsize_row = 7.5, angle_col = 0, border_color = NA, scale = "none", cluster_rows = F,keysize = 0.2,
         cellwidth = cell_width, cellheight = cell_height, width = heatmap_width, height = heatmap_height,
         display_numbers = hormone.matrix, number_color = "black")
dev.off()

### mRNA processing
cell_height <- heatmap_height/nrow(RNA.expr)
RNA.myrange <- max(abs(RNA.expr));
RNA.breaks = seq(-1.5, 1.5, length.out = 30)
RNA.heatmap <- RNA.expr[,c(1:6,13:15)]
RNA.heatmap <- RNA.heatmap[order(rowMeans(RNA.heatmap[, 1:6])), ]
RNA.matrix <- matrix("", nrow = nrow(RNA.heatmap), ncol = ncol(RNA.heatmap))
pdf("./hIAPP_pathway/mRNA_processing_hIAPP.pdf", width = 8, height = 8)
pheatmap(RNA.heatmap, main = "mRNA Processing",
         breaks = RNA.breaks,legend = T,  gaps_col = 6,labels_col = labeltiapp2dimp, 
         labels_row=toupper(rownames(RNA.heatmap)),cluster_cols = F, color = colorRampPalette(c("blue","white","red"))(30), 
         fontsize_col = 10, fontsize_row = 7.5, angle_col = 0, border_color = NA, scale = "none", cluster_rows = F,keysize = 0.2,
         cellwidth = cell_width, cellheight = cell_height, width = heatmap_width, height = heatmap_height,
         display_numbers = RNA.matrix, number_color = "black")
dev.off()

### ribosome biogenesis
cell_height <- heatmap_height/nrow(ribo.expr)
ribo.myrange <- max(abs(ribo.expr));
ribo.breaks = seq(-1.5, 1.5, length.out = 30)
ribo.heatmap <- ribo.expr[,c(1:6,13:15)]
ribo.heatmap <- ribo.heatmap[order(rowMeans(ribo.heatmap[, 1:6])), ]
ribo.matrix <- matrix("", nrow = nrow(ribo.heatmap), ncol = ncol(ribo.heatmap))
pdf("./hIAPP_pathway/ribosome_biogenesis_hIAPP.pdf", width = 8, height = 8)
pheatmap(ribo.heatmap, main = "Ribosome Biogenesis",
         breaks = ribo.breaks,legend = T,  gaps_col = 6,labels_col = labeltiapp2dimp, 
         labels_row=toupper(rownames(ribo.heatmap)),cluster_cols = F, color = colorRampPalette(c("blue","white","red"))(30), 
         fontsize_col = 10, fontsize_row = 7.5, angle_col = 0, border_color = NA, scale = "none", cluster_rows = F,keysize = 0.2,
         cellwidth = cell_width, cellheight = cell_height, width = heatmap_width, height = heatmap_height,
         display_numbers = ribo.matrix, number_color = "black")
dev.off()


