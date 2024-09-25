dds_new <- DESeq(dds)
counts <- assay(dds_new)
counts_log <- log(counts)
counts_log[is.infinite(counts_log)] <- 0
bulk.eset <- Biobase::ExpressionSet(assayData = counts_log)
sc <- read.csv("/home/kaliki_sci/Hamsters_transcriptome/Sc/sc.csv")
res <- BisqueRNA::MarkerBasedDecomposition(bulk.eset, sc, weighted=T)
marker.based.estimates <- res$bulk.props
write.csv(marker.based.estimates, "/home/kaliki_sci/Hamsters_transcriptome/Sc/PC_1_for_cell_types.csv")

#перевести в вид "тип_клеток - PC_1 - вид_червя - образец"
cell <- read.csv("/home/kaliki_sci/Hamsters_transcriptome/Sc/PC_1_for_cell_types.csv")
cell$cell_type <- as.factor(cell$cell_type)
cell$group <- as.factor(cell$group)
cell$sample <- as.factor(cell$sample)
ggplot(cell, aes(x=cell_type, y=PC_1, fill=group)) +
     geom_boxplot(position=position_dodge(0.8))+
     geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.2, 
                  position=position_dodge(0.8))
