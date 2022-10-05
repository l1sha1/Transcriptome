library(tidyverse)
library(DESeq2)
library(pheatmap)
library(org.Mm.eg.db)
library(clusterProfiler)
library(viridis)

norm_counts <- as.data.frame(counts(dds, normalized=TRUE)) # DESeq2-normalized counts: Median of ratios method https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

res <- rownames_to_column(norm_counts)

mm <- org.Mm.eg.db
entrez <- AnnotationDbi::select(mm, 
                                keys = res$rowname,
                                columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL")
entrez <- distinct(entrez, SYMBOL, .keep_all = TRUE) #removes duplicate entries if there was multiple mapping
res$entrez <- entrez$ENTREZID
write.csv(res, "/home/kaliki_sci/Hamsters_transcriptome/worms_time/entrezid_results.csv", row.names = FALSE)

dds_wt <- DESeq(dds) # test = "Wald"
OF_control <- results(dds_wt, contrast=list(c("worms_OF_vs_control")), alpha = 0.05)
write.csv(OF_control, "/home/kaliki_sci/Hamsters_transcriptome/worms_time/OF_control_Wald.csv", row.names = TRUE)

# В таблице entrezid_results.csv оставить только столбцы OF и control, добавить столбец log2FoldChanges из OF_control_Wald.csv
OF_norm_count <- read.csv("/home/kaliki_sci/Hamsters_transcriptome/worms_time/entrezid_results_OF.csv")
OF_norm_count <- filter(OF_norm_count, entrez != "NA")
geneList_OF <- OF_norm_count$log2FoldChange
names(geneList_OF) = as.character(OF_norm_count$entrez) 
gene_OF <- dplyr::filter(OF_norm_count, abs(log2FoldChange) > 1)

kk <- enrichKEGG(gene         = gene_OF,
                  organism     = 'mmu',
                  pvalueCutoff = 0.05)
kegg_df_OF <- as.data.frame(kk@result)
kegg_df_OF <- distinct(kegg_df_OF, geneID, .keep_all = TRUE)
kegg_df_20_OF <- kegg_df_OF[1:20,]
des_OF <- kegg_df_20_OF$Description
des_OF <- rev(des_OF)
kegg_df_20_OF$Description <- factor(kegg_df_20_OF$Description, levels = kegg_df_20_OF$Description[order(-(kegg_df_20_OF$pvalue))])
png(filename="/home/kaliki_sci/Hamsters_transcriptome/worms_time/OF_kegg_plot_l2fc1_BP.png", width=8, height=6, units="in", res=300)
ggplot(kegg_df_20_OF, aes(Description, Count, fill=p.adjust)) + 
     geom_col() +
     coord_flip() +
     scale_x_discrete(labels = str_wrap(des_OF, width = 60)) +
     scale_fill_viridis() +
     ylab("Number Enriched") +
     xlab("GO Term") +
     theme(
         legend.position = "right",
         text = element_text(size = 10),
         axis.title.x = element_text(size = 10), 
         axis.title.y = element_text(size = 10))
dev.off()
