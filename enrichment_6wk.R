library(enrichR)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(dplyr)

setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab/T2D")
hIAPP.6wk <- read.csv("./DEG_wk6/hIAPP_results.csv")
rIAPP.6wk <- read.csv("./DEG_wk6/rIAPP_results.csv")
hIAPP.9wk <- read.csv("./DEG_wk9/hIAPP_results.csv")
rIAPP.9wk <- read.csv("./DEG_wk9/rIAPP_results.csv")

hIAPP.6wk.up <- hIAPP.6wk[hIAPP.6wk$log2FoldChange>0 & hIAPP.6wk$padj < 0.05,]
hIAPP.6wk.down <- hIAPP.6wk[hIAPP.6wk$log2FoldChange<0 & hIAPP.6wk$padj < 0.05,]
rIAPP.6wk.up <- rIAPP.6wk[rIAPP.6wk$log2FoldChange>0 & rIAPP.6wk$padj < 0.05,]
rIAPP.6wk.down <- rIAPP.6wk[rIAPP.6wk$log2FoldChange<0 & rIAPP.6wk$padj < 0.05,]
hIAPP.9wk.up <- hIAPP.9wk[hIAPP.9wk$log2FoldChange>0 & hIAPP.9wk$padj < 0.05,]
hIAPP.9wk.down <- hIAPP.9wk[hIAPP.9wk$log2FoldChange<0 & hIAPP.9wk$padj < 0.05,]
rIAPP.9wk.up <- rIAPP.9wk[rIAPP.9wk$log2FoldChange>0 & rIAPP.9wk$padj < 0.05,]
rIAPP.9wk.down <- rIAPP.9wk[rIAPP.9wk$log2FoldChange<0 & rIAPP.9wk$padj < 0.05,]

hIAPP.6wk.up.list <- hIAPP.6wk.up$X
hIAPP.6wk.down.list <- hIAPP.6wk.down$X
rIAPP.6wk.up.list <- rIAPP.6wk.up$X
rIAPP.6wk.down.list <- rIAPP.6wk.down$X
hIAPP.9wk.up.list <- hIAPP.9wk.up$X
hIAPP.9wk.down.list <- hIAPP.9wk.down$X
rIAPP.9wk.up.list <- rIAPP.9wk.up$X
rIAPP.9wk.down.list <- rIAPP.9wk.down$X

deg.list <- list(hIAPP.6wk.up.list, rIAPP.6wk.up.list, hIAPP.9wk.up.list, rIAPP.9wk.up.list,
                 hIAPP.6wk.down.list, rIAPP.6wk.down.list, hIAPP.9wk.down.list, rIAPP.9wk.down.list)
names(deg.list) <- c("Up hIAPP vs. WT 6wk","Up rIAPP vs. WT 6wk","Up hIAPP vs. WT 9wk","Up rIAPP vs. WT 9wk",
                     "Down hIAPP vs. WT 6wk","Down rIAPP vs. WT 6wk","Down hIAPP vs. WT 9wk","Down rIAPP vs. WT 9wk")

enrich.comp <- compareCluster(geneClusters = deg.list, fun = "enrichGO", ont = "BP", 
                              OrgDb="org.Mm.eg.db", keyType = "SYMBOL", pvalueCutoff=0.05)
enrichplot::dotplot(enrich.comp, showCategory=8) + theme(axis.text.x = element_text(angle = 35, hjust = 1))

cluster_colors <- c("#F4B184", "#E6A0AD", "#99CCFF", "#A6A864", "#F4B184", "#E6A0AD", "#99CCFF", "#A6A864")
plot <- ggplot(enrich.comp, aes(x = Cluster, y = Description, color = Cluster, size = -log10(p.adjust))) +
  geom_point() +
  scale_color_manual(values = cluster_colors) +
  scale_size_continuous(breaks = c(1.3, 5, 10, 20, 30), 
                        labels = c(1.3, 5, 10, 20, 30), 
                        name = "-log10(FDR)",
                        limits = c(0, 40)) + 
  # labs(x = "Cluster", y = NULL) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        # axis.title.y = element_text(size = 14),
        plot.background = element_rect(color = "white", fill = NA),
        panel.border = element_rect(color = "black", fill = NA, size = 0.25)) +
  guides(color = FALSE)

pathway_list <- c("extracellular matrix organization", "cell chemotaxis",
                  "cell-substrate adhesion", "wound healing", "hormone secretion",
                  "mRNA processing", "ribosome biogenesis", "sterol biosynthetic process")
plot$data <- plot$data[plot$data$Description %in% pathway_list,]
levels(plot$data$Description)[1] <- "cholesterol biosynthetic process"

levels(plot$data$Cluster) <- c("Up hIAPP vs. WT 6wk","Up rIAPP vs. WT 6wk","Up hIAPP vs. WT 9wk","Up rIAPP vs. WT 9wk",
                               "Down hIAPP vs. WT 6wk","Down rIAPP vs. WT 6wk","Down hIAPP vs. WT 9wk","Down rIAPP vs. WT 9wk")
plot$data <- plot$data[plot$data$Cluster %in% c("Up hIAPP vs. WT 6wk","Up rIAPP vs. WT 6wk",
                                                "Down hIAPP vs. WT 6wk","Down rIAPP vs. WT 6wk"),]

plot

ggsave("./Enrichment/enrichment_6wk.tiff", height = 6, width = 6, units = "in")
