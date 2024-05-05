library(tximport)
library(DESeq2)

## extract data
setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab")
load("Common_Mousegenecode2gene.rda")

setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab/T2D")
hIAPP_files<-file.path("./hIAPP_raw", list.files("./hIAPP_raw"),"quant.sf")
hIAPP_files <- hIAPP_files[grep("^./hIAPP_raw/wk9", hIAPP_files)]
names(hIAPP_files) <- list.files("./hIAPP_raw")[grep("^wk9", list.files("./hIAPP_raw"))]

WT_files<-file.path("./WT_raw", list.files("./WT_raw"),"quant.sf")
WT_files <- WT_files[grep("^./WT_raw/wk9", WT_files)]
names(WT_files) <- list.files("./WT_raw")[grep("^wk9", list.files("./WT_raw"))]

files <- c(hIAPP_files, WT_files)
txi <- tximport(files, type="salmon", tx2gene=tx2gene,ignoreAfterBar = T,ignoreTxVersion = T)
head(txi$counts) # Just checking

setwd("./DEG_wk9")
write.csv(txi$counts, "hIAPP_raw_count.csv")
files_names <- names(files)

conditions <-as.data.frame(cbind(sample = c("wk9_SRX8128098", "wk9_SRX8128099", "wk9_SRX8128100",
                                            "wk9_SRX8128107", "wk9_SRX8128108", "wk9_SRX8128109"),
                                 treatment = c(rep("treatment", 3), rep("control", 3)),
                                 rownames = NULL))
write.csv(conditions, file = "hIAPP_condition.csv", row.names = FALSE)
coldata<-read.csv("hIAPP_condition.csv")
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
write.csv(res, file="hIAPP_results.csv")
