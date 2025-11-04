## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("iClusterPlus")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("genefilter")

library(magrittr)
library(reshape2)
library(iClusterPlus)
library(data.table)
library(dplyr)
library(plyr)
library(genefilter)
set.seed(123) # for reproducibility, remove for normal use

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
folds <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10')
for (fold_ in folds) {
densities <- c('', '_11', '_9', '_7', '_5', '_3')
for (density_ in densities) {

# var micro
micro_id <- paste('meta_g_taxonomic_profiles', density_, sep="")

# training
stacked_clr_all_train <- fread(
  paste('../../data/simulations/ihmp/mofa-reformatted/fold-', fold_,'-subset-train-meta_g_taxonomic_profiles', density_, "-test.tsv.gz", sep=""),
  sep="\t")
  
# seperate by view
x_metabolomics_train <- stacked_clr_all_train[stacked_clr_all_train$view == 'HMP2_metabolomics',]
x_metabolomics_train <- subset(x_metabolomics_train, select=c( "sample", "feature", "value"))
x_metabolomics_train <- xtabs(value~.,x_metabolomics_train)
x_meta_t_train <- stacked_clr_all_train[stacked_clr_all_train$view == 'meta_t_ecs',]
x_meta_t_train <- subset(x_meta_t_train, select=c( "sample", "feature", "value"))
x_meta_t_train <- xtabs(value~.,x_meta_t_train)
x_viromes_train <- stacked_clr_all_train[stacked_clr_all_train$view == 'virome_virmap_analysis',]
x_viromes_train <- subset(x_viromes_train, select=c( "sample", "feature", "value"))
x_viromes_train <- xtabs(value~.,x_viromes_train)
x_meta_g_train <- stacked_clr_all_train[stacked_clr_all_train$view == micro_id,]
x_meta_g_train <- subset(x_meta_g_train, select=c( "sample", "feature", "value"))
x_meta_g_train <- xtabs(value~.,x_meta_g_train)

# testing
stacked_clr_all_test <- fread(
  paste('../../data/simulations/ihmp/mofa-reformatted/fold-', fold_,'-subset-test-meta_g_taxonomic_profiles', density_, "-test.tsv.gz", sep=""),
  sep="\t")

# seperate by view
x_metabolomics_test <- stacked_clr_all_test[stacked_clr_all_test$view == 'HMP2_metabolomics',]
x_metabolomics_test <- subset(x_metabolomics_test, select=c( "sample", "feature", "value"))
x_metabolomics_test <- xtabs(value~.,x_metabolomics_test)
x_meta_t_test <- stacked_clr_all_test[stacked_clr_all_test$view == 'meta_t_ecs',]
x_meta_t_test <- subset(x_meta_t_test, select=c( "sample", "feature", "value"))
x_meta_t_test <- xtabs(value~.,x_meta_t_test)
x_viromes_test <- stacked_clr_all_test[stacked_clr_all_test$view == 'virome_virmap_analysis',]
x_viromes_test <- subset(x_viromes_test, select=c( "sample", "feature", "value"))
x_viromes_test <- xtabs(value~.,x_viromes_test)
x_meta_g_test <- stacked_clr_all_test[stacked_clr_all_test$view == micro_id,]
x_meta_g_test <- subset(x_meta_g_test, select=c( "sample", "feature", "value"))
x_meta_g_test <- xtabs(value~.,x_meta_g_test)

# compile these into a single X object
x_metabolomics_test <- as.data.frame.matrix(x_metabolomics_test)
x_metabolomics_train <- as.data.frame.matrix(x_metabolomics_train)
x_metabolomics <- merge(x_metabolomics_test, x_metabolomics_train, all = TRUE)
x_metabolomics[is.na(x_metabolomics)] <- 0
rownames(x_metabolomics) <-  c(rownames(x_metabolomics_test), rownames(x_metabolomics_train))

x_meta_t_test <- as.data.frame.matrix(x_meta_t_test)
x_meta_t_train <- as.data.frame.matrix(x_meta_t_train)
x_meta_t <- merge(x_meta_t_test ,x_meta_t_train, all = TRUE)
x_meta_t[is.na(x_meta_t)] <- 0
rownames(x_meta_t) <-  c(rownames(x_meta_t_test), rownames(x_meta_t_train))

x_viromes_test <- as.data.frame.matrix(x_viromes_test)
x_viromes_train <- as.data.frame.matrix(x_viromes_train)
x_viromes <- merge(x_viromes_test, x_viromes_train, all = TRUE)
x_viromes[is.na(x_viromes)] <- 0
rownames(x_viromes) <-  c(rownames(x_viromes_test), rownames(x_viromes_train))

x_meta_g_test <- as.data.frame.matrix(x_meta_g_test)
x_meta_g_train <- as.data.frame.matrix(x_meta_g_train)
x_meta_g <- merge(x_meta_g_test, x_meta_g_train, all = TRUE)
x_meta_g[is.na(x_meta_g)] <- 0
rownames(x_meta_g) <-  c(rownames(x_meta_g_test), rownames(x_meta_g_train))

X <- list(metabolomics = as.data.frame.matrix(x_metabolomics),
          metatranscriptomics = as.data.frame.matrix(x_meta_t),
          viromes = as.data.frame.matrix(x_viromes),
          metagenomics = as.data.frame.matrix(x_meta_g))

# make sure names are shared
shared_names <- unique(Reduce(intersect,  lapply(X, rownames)))

metabolomics.data <- as.data.frame.matrix(x_metabolomics[rownames(x_metabolomics) %in% shared_names, ])
metatranscriptomics.data <- as.data.frame.matrix(x_meta_t[rownames(x_meta_t) %in% shared_names, ])
viromes.data <- as.data.frame.matrix(x_viromes[rownames(x_viromes) %in% shared_names, ])
metagenomics.data <- as.data.frame.matrix(x_meta_g[rownames(x_meta_g) %in% shared_names, ])

ordered_names <- order(rownames(metabolomics.data))
ordered_names_save <- rownames(metabolomics.data)[ordered_names]
metabolomics.data <-metabolomics.data[ordered_names,]
metatranscriptomics.data <-metatranscriptomics.data[ordered_names,]
viromes.data <-viromes.data[ordered_names,]
metagenomics.data <-metagenomics.data[ordered_names,]

metabolomics.data <- data.matrix(metabolomics.data)
metatranscriptomics.data <- data.matrix(metatranscriptomics.data)
viromes.data <- data.matrix(viromes.data)
metagenomics.data <- data.matrix(metagenomics.data)

# have to reduce data size for iCuster (following the tutorial)
# https://bioconductor.org/packages/release/bioc/vignettes/iClusterPlus/inst/doc/iManual.pdf
metabolomics.data.lv <- t(varFilter(t(metabolomics.data), var.cutoff=0.95))
metatranscriptomics.data.lv <- t(varFilter(t(metatranscriptomics.data), var.cutoff=0.95))

# run iCluster
fit.single <- iClusterPlus(dt1=metabolomics.data.lv, dt2=metatranscriptomics.data.lv, dt3=viromes.data, dt4=metagenomics.data,
                           type=c("gaussian","gaussian","gaussian","gaussian"),
                           lambda=c(0.2,0.2,0.2,0.2),
                           K=3, maxiter=10)

loadings.out <- as.data.frame.matrix(fit.single$meanZ)
rownames(loadings.out) <- ordered_names_save
write.csv(loadings.out, paste("../../results/simulations/ihmp/icluster/", fold_ ,".factors.model.", density_, ".csv", sep=""))

}
}
