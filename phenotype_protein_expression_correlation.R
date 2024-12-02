#Clear R environment
rm(list = ls())

#Import differential expression data
reserve <- read.csv(file.choose(), header=TRUE)

#Proteins to fixed row identifiers
rownames(reserve) = reserve$protein
reserve$protein = NULL

#Correlation of differential protein expression between phenotypes 
library(stats)
library(RVAideMemoire)

cor.test(reserve$logfc_reactive_uninflamed, reserve$logfc_hyper_hypo, method = "spearman")
spearman.ci(reserve$logfc_reactive_uninflamed, reserve$logfc_hyper_hypo, nrep = 10000, conf.level = 0.95)
