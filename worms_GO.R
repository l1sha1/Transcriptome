OF_control <- results(dds_wt, contrast=list(c("worms_OF_vs_control")), alpha = 0.05)
write.csv(OF_control, "/home/kaliki_sci/Hamsters_transcriptome/worms_time/OF_control_Wald.csv", row.names = TRUE)
OF_norm_count <- read.csv("/home/kaliki_sci/Hamsters_transcriptome/worms_time/entrezid_results_OF.csv")
OF_norm_count <- filter(OF_norm_count, entrez != "NA")
geneList <- OF_norm_count$log2FoldChange
names(geneList) = as.character(OF_norm_count$entrez) 
gene <- dplyr::filter(OF_norm_count, abs(log2FoldChange) > 2)

write.csv(gene, "/home/kaliki_sci/Hamsters_transcriptome/worms_time/OF_sig_genes.csv", row.names = FALSE)
gene <- gene$entrez

ego <- enrichGO(gene          = gene,
universe      = names(geneList),
OrgDb         = org.Mm.eg.db,
ont           = "BP",       
pAdjustMethod = "BH",
pvalueCutoff  = 0.05,
qvalueCutoff  = 0.1,
readable      = TRUE)
go_df <- as.data.frame(ego@result)

go_df <- distinct(go_df, geneID, .keep_all = TRUE)
go_df_20 <- go_df[1:20,]
des <- go_df_20$Description
des <- rev(des)
go_df_20$Description <- factor(go_df_20$Description, levels = go_df_20$Description[order(-(go_df_20$pvalue))])
png(filename="/home/kaliki_sci/Hamsters_transcriptome/worms_time/OF_GO_plot_l2fc2_BP.png", width=8, height=6, units="in", res=300)
ggplot(go_df_20, aes(Description, Count, fill=p.adjust)) + 
geom_col() +
coord_flip() +
scale_x_discrete(labels = str_wrap(des, width = 60)) +
scale_fill_viridis() +
ylab("Number Enriched") +
xlab("GO Term") +
theme(
legend.position = "right",
text = element_text(size = 10),
axis.title.x = element_text(size = 10), 
axis.title.y = element_text(size = 10))
> dev.off()
