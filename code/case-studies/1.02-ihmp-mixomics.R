## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(magrittr)
library(reshape2)
library(mixOmics)
library(data.table)
library(dplyr)
library(plyr)
set.seed(123) # for reproducibility, remove for normal use

# train
stacked_clr_all_train <- fread(
  '../../data/case-studies/ihmp/mofa-tables/subset-train.tsv.gz',
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
x_meta_g_train <- stacked_clr_all_train[stacked_clr_all_train$view == 'meta_g_taxonomic_profiles',]
x_meta_g_train <- subset(x_meta_g_train, select=c( "sample", "feature", "value"))
x_meta_g_train <- xtabs(value~.,x_meta_g_train)
x_prot_train <- stacked_clr_all_train[stacked_clr_all_train$view == 'HMP2_proteomics_ecs',]
x_prot_train <- subset(x_prot_train, select=c( "sample", "feature", "value"))
x_prot_train <- xtabs(value~.,x_prot_train)

# testing
stacked_clr_all_test <- fread(
  '../../data/case-studies/ihmp/mofa-tables/subset-test.tsv.gz',
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
x_meta_g_test <- stacked_clr_all_test[stacked_clr_all_test$view == 'meta_g_taxonomic_profiles',]
x_meta_g_test <- subset(x_meta_g_test, select=c( "sample", "feature", "value"))
x_meta_g_test <- xtabs(value~.,x_meta_g_test)
x_prot_test <- stacked_clr_all_test[stacked_clr_all_test$view == 'HMP2_proteomics_ecs',]
x_prot_test <- subset(x_prot_test, select=c( "sample", "feature", "value"))
x_prot_test <- xtabs(value~.,x_prot_test)

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


x_prot_test <- as.data.frame.matrix(x_prot_test)
x_prot_train <- as.data.frame.matrix(x_prot_train)
x_prot <- merge(x_prot_test, x_prot_train, all = TRUE)
x_prot[is.na(x_prot)] <- 0
rownames(x_prot) <-  c(rownames(x_prot_test), rownames(x_prot_train))


X <- list(metabolomics = as.data.frame.matrix(x_metabolomics),
          proteomics = as.data.frame.matrix(x_prot),
          metatranscriptomics = as.data.frame.matrix(x_meta_t),
          viromes = as.data.frame.matrix(x_viromes),
          metagenomics = as.data.frame.matrix(x_meta_g))

# make sure names are shared
shared_names <- unique(Reduce(intersect,  lapply(X, rownames)))
X <- list(metabolomics = as.data.frame.matrix(x_metabolomics[rownames(x_metabolomics) %in% shared_names, ]),
          proteomics = as.data.frame.matrix(x_prot[rownames(x_prot) %in% shared_names, ]),
          metatranscriptomics = as.data.frame.matrix(x_meta_t[rownames(x_meta_t) %in% shared_names, ]),
          viromes = as.data.frame.matrix(x_viromes[rownames(x_viromes) %in% shared_names, ]),
          metagenomics = as.data.frame.matrix(x_meta_g[rownames(x_meta_g) %in% shared_names, ]))



## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DIABLO in mixOmics with a value of 0.5 in the design matrix 
# (cite: https://diabetesjournals.org/care/article/41/10/2178/36703/Intestinal-Metaproteomics-Reveals-Host-Microbiota)
design = matrix(0.5, 
                ncol = length(X),
                nrow = length(X),
                dimnames = list(names(X),
                                names(X)))
diag(design) = 0 # set diagonal to 0s
model.spls <- block.spls(X, indY=5, ncomp = 3, design = 1, scale = TRUE)
factors <- ldply(Map(cbind, model.spls$variates, subject = lapply(model.spls$variates, rownames)))
write.csv(factors,  '../../data/case-studies/ihmp/other-tools/mixomics-loadings.csv')
