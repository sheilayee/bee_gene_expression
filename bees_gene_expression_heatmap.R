library(ggplot2)
library(pheatmap)
library(dplyr)

setwd('/Users/Sheila/Documents/BF527/project/')
pwd=getwd()

# read in normalized counts
norm_counts <- read.csv('/Users/Sheila/Documents/BF527/project/normalized_counts.csv', row.names = 1)

# read in ordered output of DESeq2
deseq_res <- read.csv('/Users/Sheila/Documents/BF527/project/normalized_ordered_counts.csv', row.names = 1)

# load metadata
metadata <- read.csv('/Users/Sheila/Documents/BF527/project/Design.txt', header = T, sep =  "\t")
metadata$Gender <- factor(metadata$Gender, levels=c('Male', 'Female'), labels=c("Male", "Female"))
str(metadata$Gender)
conditions <- data.frame(metadata[c("Gender", "Caste", "DevStage")])
row.names(conditions) <- metadata$X
conditions

# add in column indicating if p-adjusted is less than 0.05
# deseq_res$sig <- ifelse(deseq_res$padj <= 0.05, "yes", "no")
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
# ggplot(deseq_res, aes(x = log2FoldChange, y = -log10(padj), color = sig)) + 
#   geom_point()
p <- ggplot(deseq_res, aes(x = log2FoldChange, y = -log10(padj), color = significant, label= significant)) + 
  geom_point(size=0.5) +
  theme_bw() + 
  scale_color_manual(values=c("gray70", "steelblue1", "firebrick")) +
  ggtitle("Volcano plot comparing \n male and female samples") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.text=element_text(size=10)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "goldenrod1", size=0.5) + 
  geom_vline(xintercept=-log2(2), linetype="dashed", color = "goldenrod1", size=0.5) +
  geom_vline(xintercept=log2(2), linetype="dashed", color = "goldenrod1", size=0.5)
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
sig_norm_counts <- sig_norm_counts[,2:28]
row.names(sig_norm_counts) <- sig_gene_names

# another way to do it
sig_gene_names <- rownames(sig_genes)
sig_norm_counts2 = norm_counts[which(rownames(norm_counts) %in% sig_gene_names), ]

# plot heatmap
ann_color = list("Gender"=c("Male"="lightskyblue","Female"="maroon1"),
                 "Caste"=c("male"="skyblue3","queen"="olivedrab3", "worker"="darksalmon"), 
                 "DevStage"=c("Adult"="mediumslateblue", "Larva"="lightgoldenrod1", "Pupa"="plum2"))


pheatmap(log2(sig_norm_counts2+1), scale='row', show_rownames = F, treeheight_row = 0, 
         annotation_col = conditions, 
         annotation_colors = ann_color)
savePlot(filename, type) 


