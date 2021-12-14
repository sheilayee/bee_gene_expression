# install necessary packages
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")

# BiocManager::install("DESeq2", force = TRUE)
# install.packages("tidyverse")
# install.packages("pheatmap")
# install.packages("dplyr")

library(pheatmap)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(tibble)
library(ggplot2)

setwd('/Users/Sheila/Documents/BF527/project/')
pwd=getwd()

# load count matrix
# this txt file has no genes with less than 50 counts
# subset the file to only consider samples that are labeled either reproductive or non-reproductive
counts <- (read.csv('/Users/Sheila/Documents/BF527/project/Counts.txt', sep =  "\t", header = T, row.names = 1)[,-c(1:15,25:26)])

# load metadata
# subset the file to only consider samples that are labeled either reproductive or non-reproductive
metadata <- (read.csv('/Users/Sheila/Documents/BF527/project/Design.txt', header = T, sep =  "\t")[-c(1:15,25:26),])

# change transcript column name into row headers
# row.names(counts) <- counts$Transcript
# # counts <- round(apply(counts,2,as.numeric))
# counts$Transcript <- NULL

# remove missing counts data
counts <- na.omit(counts)

# change Gender column into factor
# metadata <- metadata[16:27,]
metadata$Reproductive <- factor(metadata$Reproductive, levels=c('N', 'R'), labels=c("Non_reproductive", "Reproductive"))
metadata$Caste <- factor(metadata$Caste, levels=c('worker', 'queen'), labels=c("worker", "queen"))
str(metadata$Reproductive)
str(metadata$Caste)
metadata$Reproductive

# create data object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~Caste + Reproductive)

# run DESeq2 (this step includes normalization step)
# we have to take into account the number of reads mapped to different sized libraries
dds <- DESeq(dds)

# extract normalized read counts
norm_counts <- counts(dds, normalized = T)
write.csv(norm_counts, "normalized_counts_repro_non-repro.csv")

# look at DESeq results
# only consider genes with adjusted p-value of less than 0.05 (that is the 
# threshold it must meet to be a DEG)
res <- results(dds, alpha = 0.05)

# look at what the conditions being compared are  
resultsNames(dds)

# summary of differential gene expression
# tells you how many genes were measured, the corrected p-value level, 
# LFC is log fold change, and the outliers
summary(res)

# sort summary list by smallest to largest adjusted p-value
res_ordered <- res[order(res$padj),]
write.csv(res_ordered, "normalized_ordered_counts_repro_non-repro.csv")
head(res_ordered)
tail(res_ordered)


