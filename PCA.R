library(DESeq2)
library(ggforce)
library(cowplot)
library(readxl)

setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab/T2D/")

data.wk6 <- read.csv("./DEG_wk6/hIAPP_raw_count.csv")
data.wk9 <- read.csv("./DEG_wk9/hIAPP_raw_count.csv")
data <- merge(data.wk6, data.wk9,by.x = "X",by.y = "X")
rownames(data) <- data$X
data <- data[,-1]
data <- as.matrix(data)
data <- round(data)

target.data <- as.data.frame(c(rep("hIAPP_6wk",6),rep("WT_6wk",6),rep("hIAPP_9wk",3), rep("WT_9wk",3)))
colnames(target.data) <- "Group"

dds <- DESeqDataSetFromMatrix(data,colData = target.data, design = ~ Group)

vsd <- vst(dds)

sample.colors <- c("hIAPP_6wk"="#B856D7","WT_6wk" = "#89A6D2","hIAPP_9wk" = "#A0A0A4","WT_9wk" = "#D4D4D4")
sample.shape <- c("hIAPP_6wk" = 16, "WT_6wk" = 1, "hIAPP_9wk" = 17, "WT_9wk" = 2)

pcaData <- plotPCA(vsd, intgroup="Group", ntop = dim(data)[[1]], returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"),digits=2)

setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab/T2D/PCA")
tiff("PCA_1.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, color = group, shape = group)) +
  geom_point(size = 2.5) +
  geom_mark_ellipse(aes(fill = NULL, color = group)) +
  scale_colour_manual(values = sample.colors, breaks = names(sample.colors)) +
  scale_shape_manual(values = sample.shape, breaks = names(sample.shape)) +
  theme_half_open(font_size = 15) +
  xlab(paste0("Dimension 1")) +
  ylab(paste0("Dimension 2")) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.1, 0.2),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  ylim(-40, 40) +
  xlim(-40, 40) +
  theme(axis.text = element_text(size = 11))
dev.off()



