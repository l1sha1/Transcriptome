library('msigdb')
mm_msigdb <- getMsigdb('mm', 'EZID')
listCollections(mm_msigdb)
# "c1" "c3" "c2" "c8" "c6" "c7" "c4" "c5" "h" 
    # H: hallmark gene sets
    # C1: positional gene sets
    # C2: curated gene sets
    # C3: motif gene sets
    # C4: computational gene sets
    # C5: GO gene sets
    # C6: oncogenic signatures
    # C7: immunologic signatures
    
h = subsetCollection(mm_msigdb, 'h')
h_gmt <- tempfile()
toGmt(h, h_gmt)
h_1 <- read.gmt(h_gmt)
h_OF <- enricher(gene_OF, TERM2GENE=h_1)
msigdb_h_OF <- as.data.frame(h_OF@result)

write.csv(msigdb_h_OF, "/home/kaliki_sci/Hamsters_transcriptome/worms_time/MSIGDB_H_OF_l2fc1.csv")

msigdb_h_OF <- distinct(msigdb_h_OF, geneID, .keep_all = TRUE)
msigdb_h_OF_20 <- msigdb_h_OF[1:20,]
des_OF_h <- msigdb_h_OF_20$Description
des_OF_h <- rev(des_OF_h)
msigdb_h_OF_20$Description <- factor(msigdb_h_OF_20$Description, levels = msigdb_h_OF_20$Description[order(-(msigdb_h_OF_20$pvalue))])
png(filename="/home/kaliki_sci/Hamsters_transcriptome/worms_time/OF_MSIGDB_H_plot_l2fc1_BP.png", width=8, height=6, units="in", res=300)

ggplot(msigdb_h_OF_20, aes(Description, Count, fill=p.adjust)) + 
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = str_wrap(des_OF_h, width = 60)) +
  scale_fill_viridis() +
  ylab("Number Enriched") +
  xlab("GO Term") +
  theme(
    legend.position = "right",
    text = element_text(size = 10),
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10))
dev.off()
