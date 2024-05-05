library(tximport)
library(DESeq2)

## extract data
setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab")
load("Common_Mousegenecode2gene.rda")

setwd("/Users/rainyliu/Desktop/Learning/UCLA_Yang_Lab/T2D")
rIAPP_files<-file.path("./rIAPP_raw", list.files("./rIAPP_raw"),"quant.sf")
rIAPP_files <- rIAPP_files[grep("^./rIAPP_raw/wk9", rIAPP_files)]
names(rIAPP_files) <- list.files("./rIAPP_raw")[grep("^wk9", list.files("./rIAPP_raw"))]

WT_files<-file.path("./WT_raw", list.files("./WT_raw"),"quant.sf")
WT_files <- WT_files[grep("^./WT_raw/wk9", WT_files)]
names(WT_files) <- list.files("./WT_raw")[grep("^wk9", list.files("./WT_raw"))]

files <- c(rIAPP_files, WT_files)
txi <- tximport(files, type="salmon", tx2gene=tx2gene,ignoreAfterBar = T,ignoreTxVersion = T)
head(txi$counts) # Just checking

setwd("./DEG_wk9")
write.csv(txi$counts, "rIAPP_raw_count.csv")
files_names <- names(files)

conditions <-as.data.frame(cbind(sample = c("wk9_SRX8128104", "wk9_SRX8128105", "wk9_SRX8128106",
                                            "wk9_SRX8128107", "wk9_SRX8128108", "wk9_SRX8128109"),
                                 treatment = c(rep("treatment", 3), rep("control", 3)),
                                 rownames = NULL))
write.csv(conditions, file = "rIAPP_condition.csv", row.names = FALSE)
coldata<-read.csv("rIAPP_condition.csv")
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
write.csv(res, file="rIAPP_results.csv")
