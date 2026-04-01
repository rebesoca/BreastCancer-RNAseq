args <- commandArgs(trailingOnly=TRUE)
counts_file <- args[1]

library(DESeq2)

counts <- read.table(counts_file, header=TRUE, row.names=1)

counts <- counts[, -c(1:5)]

condition <- factor(c(rep("BTZ",4), rep("DOX",4), rep("CTR",4)))
coldata <- data.frame(row.names=colnames(counts), condition)

dds <- DESeqDataSetFromMatrix(countData = counts,
                             colData = coldata,
                             design = ~ condition)

dds <- DESeq(dds)

res <- results(dds)

write.csv(as.data.frame(res), "DESeq2_results.csv")
