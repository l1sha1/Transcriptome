#!/usr/bin/env Rscript

library("GenomicFeatures")
library("DESeq2")
library("dplyr")
library("ggplot2")
library("GenomicAlignments")
library("BiocParallel")
library("RColorBrewer")
library("pheatmap")


gtffile <- file.path("/beegfs/scratch/ws/ws1/lishaiea-work/genomic.gff") ##file_name
txdb <- makeTxDbFromGFF(gtffile, circ_seqs = character())
ebg <- exonsBy(txdb, by="gene")

csvfile <- file.path("/beegfs/scratch/ws/ws1/lishaiea-work/sample_table.csv")
sampleTable <- read.csv(csvfile, header = T)
sampleTable$Condition <- as.factor(sampleTable$Condition)
sampleTable$Condition <- relevel(sampleTable$Condition, ref = 'control')


print(sampleTable)
filenames <- file.path("/beegfs/scratch/ws/ws1/lishaiea-work/H69", paste0(sampleTable$Sample))
file.exists(filenames)
print('files exist')
se <- summarizeOverlaps(features=ebg, reads=filenames, mode="Union", singleEnd=FALSE,  ignore.strand=TRUE, fragments=TRUE)

colData(se) <- DataFrame(sampleTable)

#how many maped reads to sample?
round( colSums(assay(se)) / 1e6, 1 )
print('mapped reads to sample')
dds <- DESeqDataSet(se, design = ~ Condition)

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 24, ]
nrow(dds)


#rld <- rlog(dds, blind = FALSE)
#vsd <- vst(dds, blind = FALSE)
#head(assay(rld), 3)
#head(assay(vsd), 3)
print("passed, count dds")
dds <- estimateSizeFactors(dds)
head(dds)
dds_new <- DESeq(dds)
res <- results(dds_new)
summary(res)
head(res)
write.csv(res[order(res$padj), ], file = 'deseq_result_star_SA.csv')
save(dds, file="dds_h69.RData")
save(se, file="countMat_h69.RData")

