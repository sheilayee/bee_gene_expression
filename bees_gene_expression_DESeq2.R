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
counts <- read.csv('/Users/Sheila/Documents/BF527/project/Counts.txt', sep =  "\t", header = T, row.names = 1)

# load metadata
metadata <- read.csv('/Users/Sheila/Documents/BF527/project/Design.txt', header = T, sep =  "\t")

# change transcript column name into row headers
# row.names(counts) <- counts$Transcript
# # counts <- round(apply(counts,2,as.numeric))
# counts$Transcript <- NULL

# count how many transcripts there are
nrow(counts)

# remove missing counts data
counts <- na.omit(counts)

# change Gender column into factor
metadata$Gender <- factor(metadata$Gender, levels=c('Male', 'Female'), labels=c("Male", "Female"))
# metadata$Caste <- factor(metadata$Caste, levels=c('male', 'worker', "queen"), labels=c("male", "worker", "queen"))
metadata$DevStage <- factor(metadata$DevStage, levels=c('Larva', 'Pupa', 'Adult'), labels=c('Larva', 'Pupa', 'Adult'))
str(metadata$Gender)
str(metadata$Caste)
str(metadata$DevStage)
metadata$Caste

# change male caste factor (to allow for linear regression in DESeq2)
metadata$Caste <- as.character(metadata$Caste)
metadata$Caste[metadata$Caste == "male"] <- "m"
metadata$Caste <- factor(metadata$Caste, levels=c('m', 'worker', "queen"), labels=c("m", "worker", "queen"))

# create data object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~DevStage + Gender)

# run DESeq2 (this step includes normalization step)
# we have to take into account the number of reads mapped to different sized libraries
dds <- DESeq(dds)

# extract normalized read counts
norm_counts <- counts(dds, normalized = T)
write.csv(norm_counts, "normalized_counts.csv")

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
write.csv(res_ordered, "normalized_ordered_counts.csv")
head(res_ordered)
tail(res_ordered)


