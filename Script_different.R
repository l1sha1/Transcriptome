# тепловая карта для всех генов
heatmap.2(counts_norm, margins = c(6,15), scale="row", trace="none", Colv = NA, Rowv = NA, dendrogram = "none",col = redgreen(75), density.info="none", cexRow=0.6)
# тепловая карта для 30 генов с наибольшим изменением Z-score
topVarGenes <- head(order(rowVars(counts_norm), decreasing=TRUE), 30)

pheatmap(CC, scale = "row", fontsize_row = 5, fontsize_col = 7, cluster_rows = TRUE, cluster_cols = FALSE, legend = TRUE)


heatmap.2(counts_norm[topVarGenes,], margins = c(6,15), scale="row", trace="none", Colv = NA, Rowv = NA, dendrogram = "none",col = redgreen(75), density.info="none", cexRow=1)
# dotplot gene-group
data_bp %>% 
     filter(data_bp$norm_gene_expression > 0) %>%
     ggplot(aes(x=group, y = cluster, size = norm_gene_expression)) + 
     theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
     geom_point() 
