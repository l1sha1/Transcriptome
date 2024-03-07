library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)
CS_l2fc1_padj0_05

CS <- read.csv("/home/kaliki_sci/Hamsters_transcriptome/worms_time/CS_control_Wald.csv")
CS <- as.data.frame(CS)
CS$significant <- ifelse(CS$padj < .05 & abs(CS$log2FoldChange) > 1 , "Significant", NA)

ggplot(CS, aes(log2FoldChange, -log10(padj), colour=significant))+ geom_point(size=1) + scale_y_continuous(limits=c(0, 15), oob=squish)+  scale_x_continuous(limits=c(-9, 10), oob=squish)

ggplot(CS, aes(baseMean, log2FoldChange, colour=significant)) + geom_point(size=1) + scale_y_continuous(limits=c(-10, 10), oob=squish) + scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=1) + labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()



# esp_control <- results(dds, contrast=list(c("Condition_esp_vs_control")))
# esp <- as.data.frame(esp_control)
esp$significant <- "NO"
esp$significant[esp$log2FoldChange > 0.6 & esp$padj < 0.05] <- "UP"
esp$significant[esp$log2FoldChange < -0.6 & esp$padj < 0.05] <- "DOWN"
ggplot(data = esp, aes(x = log2FoldChange, y = -log10(padj), col = significant)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_point(size = 1) + scale_y_continuous(limits=c(0, 40), oob=squish) +  scale_x_continuous(limits=c(-5, 5), oob=squish) +
    scale_color_manual(values = c("#FF0000", "grey", "#00CC00"),
                       labels = c("Downregulated", "Not significant", "Upregulated")) 
