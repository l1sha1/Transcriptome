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

write.csv(gene, "/home/kaliki_sci/Hamsters_transcriptome/worms_time/OF_sig_genes.csv", row.names = FALSE)
gene_OF <- gene_OF$entrez

ego_OF <- enrichGO(gene          = gene_OF,
                universe      = names(geneList_OF),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",        #BP, CC, or MF
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)
go_df_OF <- as.data.frame(ego_OF@result)
write.csv(go_df_OF, "/home/kaliki_sci/Hamsters_transcriptome/worms_time/GO_OV_l2fc1.csv")
go_df_OF <- distinct(go_df_OF, geneID, .keep_all = TRUE)
go_df_20_OF <- go_df_OF[1:20,]
des_OF <- go_df_20_OF$Description
des_OF <- rev(des_OF)
go_df_20_OF$Description <- factor(go_df_20_OF$Description, levels = go_df_20_OF$Description[order(-(go_df_20_OF$pvalue))])
png(filename="/home/kaliki_sci/Hamsters_transcriptome/worms_time/OF_GO_plot_l2fc2_BP.png", width=8, height=6, units="in", res=300)
# png(filename="/home/kaliki_sci/Hamsters_transcriptome/worms_time/OF_GO_plot_l2fc1_BP.png", width=8, height=6, units="in", res=300), если abs(log2FC)> 1
ggplot(go_df_20_OF, aes(Description, Count, fill=p.adjust)) + 
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
