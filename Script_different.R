# тепловая карта для всех генов
heatmap.2(counts_norm, margins = c(6,15), scale="row", trace="none", Colv = NA, Rowv = NA, dendrogram = "none",col = redgreen(75), density.info="none", cexRow=0.6)
# тепловая карта для 30 генов с наибольшим изменением Z-score
topVarGenes <- head(order(rowVars(counts_norm), decreasing=TRUE), 30)
heatmap.2(counts_norm[topVarGenes,], margins = c(6,15), scale="row", trace="none", Colv = NA, Rowv = NA, dendrogram = "none",col = redgreen(75), density.info="none", cexRow=1)
