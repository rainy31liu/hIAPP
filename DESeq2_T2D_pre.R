library(tximport)
library(DESeq2)

## extract data
setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab")
load("Common_Humangenecode2gene.rda")

setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab/T2D")
pre_files<-file.path("./T2D_raw/pre", list.files("./T2D_raw/pre"),"quant.sf")
names(pre_files) <- list.files("./T2D_raw/pre")

WT_files<-file.path("./T2D_raw/WT", list.files("./T2D_raw/WT"),"quant.sf")
names(WT_files) <- list.files("./T2D_raw/WT")

files <- c(pre_files, WT_files)
txi <- tximport(files, type="salmon", tx2gene=tx2gene,ignoreAfterBar = T,ignoreTxVersion = T)
head(txi$counts) # Just checking

setwd("./DEG_T2D")
write.csv(txi$counts, "pre_raw_count.csv")
files_names <- names(files)
conditions <-as.data.frame(cbind(sample = c(files_names[grepl("^pre",files_names)],
                                            files_names[grep("^WT",files_names)]),
                                 treatment = c(rep("treatment", sum(grepl("^pre",files_names))), 
                                               rep("control", sum(grepl("^WT",files_names)))),
                                 rownames = NULL))

write.csv(conditions, file = "pre_conditions.csv", row.names = FALSE)
coldata<-read.csv("pre_conditions.csv")
colnames(coldata) <- c("sample", "treatment")

#running deseq2
dds <- DESeqDataSetFromTximport(txi,colData=coldata,design = ~treatment)
keep <- rowMeans(counts(dds) >= 1) >= 0.5
dds <- dds[keep,]
dds <- DESeq(dds)
res <- as.data.frame(results(dds, contrast = c("treatment", "treatment", "control")))
summary(res)

# Order by p-value
res <- res[order(res$pvalue),]

# Saving as csv
write.csv(res, file="pre_results.csv")
