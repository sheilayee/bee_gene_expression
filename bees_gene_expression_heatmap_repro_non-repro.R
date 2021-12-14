library(ggplot2)
library(pheatmap)
library(dplyr)


setwd('/Users/Sheila/Documents/BF527/project/')
pwd=getwd()

# read in normalized counts
norm_counts <- read.csv('/Users/Sheila/Documents/BF527/project/normalized_counts_repro_non-repro.csv', row.names = 1)

# read in ordered output of DESeq2
deseq_res <- read.csv('/Users/Sheila/Documents/BF527/project/normalized_ordered_counts_repro_non-repro.csv', row.names = 1)

# load metadata
metadata <- (read.csv('/Users/Sheila/Documents/BF527/project/Design.txt', header = T, sep =  "\t")[-c(1:15,25:26),])
metadata$Reproductive <- factor(metadata$Reproductive, levels=c('N', 'R'), labels=c("Non-reproductive", "Reproductive"))
str(metadata$Reproductive)
metadata$Reproductive
conditions <- data.frame(metadata[c("Caste","Reproductive")])
row.names(conditions) <- metadata$X
conditions

# add in column indicating if p-adjusted is less than 0.05
deseq_res$significant <- ifelse(deseq_res$padj <= 0.05 & deseq_res$log2FoldChange > 0, "Sig. up-regulated", 
                                ifelse(deseq_res$padj <= 0.05 & deseq_res$log2FoldChange < 0, "Sig. down-regulated",
                                       "Not significant"))

# remove any na's if applicable
deseq_res <- na.omit(deseq_res)

# replicating built-in function in DESe12 called plotMA
# show the relationship b/w the logfoldchange of each gene and the number of 
# normalized reads. will also show if a gene is sig or not
ggplot(deseq_res, aes(x = log10(baseMean), y = log2FoldChange, color = significant)) +
  geom_point()

# volcano plot
# y-axis is indicating how many 0's is in adj p-value (higher on the y-axis, the more significant)
theme_update(plot.title = element_text(hjust = 0.5))

p <- ggplot(deseq_res, aes(x = log2FoldChange, y = -log10(padj), color = significant, label= significant)) + 
  geom_point(size=0.5) +
  theme_bw() + 
  scale_color_manual(values=c("gray70", "steelblue1", "firebrick")) +
  ggtitle("Volcano plot of samples with \n reproductive status") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text=element_text(size=10)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "goldenrod1", size=0.7) + 
  geom_vline(xintercept=-log2(2), linetype="dashed", color = "goldenrod1", size=0.7) +
  geom_vline(xintercept=log2(2), linetype="dashed", color = "goldenrod1", size=0.7)
p + guides(color = guide_legend(override.aes = list(size = 2)))



# heatmap
# look at normalized read count data of the significantly differentially expressed genes

# determine significantly expressed genes
sig_genes <- deseq_res[!is.na(deseq_res$padj) &
                         deseq_res$padj<0.05 &
                   abs(deseq_res$log2FoldChange)>=1,]

# extract the significant gene names
sig_gene_names <- rownames(sig_genes)

# merge deseq results and normalized read count using row names
sig_norm_counts <- merge(norm_counts, sig_genes, by = 0)
sig_norm_counts <- sig_norm_counts[,2:11]
row.names(sig_norm_counts) <- sig_gene_names

# another way to do it
sig_gene_names <- rownames(sig_genes)
sig_norm_counts2 = norm_counts[which(rownames(norm_counts) %in% sig_gene_names), ]


d <- dist( t(sig_norm_counts2) , method="euclidean")
sample_cor <- cor( sig_norm_counts2 )
round(sample_cor,4)

#Transform the scale from correlations
cor_distance <- -(sample_cor - 1)/2
round(cor_distance,4)

#Convert it to a distance object
d2 <- as.dist(cor_distance)
d2

# plot heatmap
ann_color = list("Caste"=c("queen"="olivedrab3", "worker"="darksalmon"), 
                 "Reproductive"=c("Non-reproductive"="lightskyblue", "Reproductive"="plum2"))


pheatmap(log2(sig_norm_counts2+1), scale='row', show_rownames = F, treeheight_row = 0, 
         annotation_col = conditions, 
         annotation_colors = ann_color)
# savePlot(filename, type) 


