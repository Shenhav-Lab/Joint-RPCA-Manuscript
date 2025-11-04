## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("mixOmics")

library(magrittr)
library(reshape2)
library(mixOmics)
library(data.table)
library(dplyr)
library(plyr)
set.seed(123) # for reproducibility, remove for normal use

setwd('/Users/bec5786/Desktop/Shenhav Lab/github/Joint-RPCA-Manuscript/code/simulations-benchmarking/')

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
folds <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10')
fsignals <- c('5', '10', '20', '30', '40', '50', '60')

for (fold_ in folds) {
  
  for (nfeats in fsignals) {
    # training
    stacked_clr_all_train <- fread(
      paste('../../data/simulations/lowrank/mofa_tables_fsignal/fold-', fold_,
            '-subset-train-fsignal', nfeats, ".tsv.gz", sep=""), sep="\t")
      
    # separate by view (# of features in each modality)
    x_f100_train <- stacked_clr_all_train[stacked_clr_all_train$view == 'f_100',]
    x_f100_train <- subset(x_f100_train, select=c("sample", "feature", "value"))
    x_f100_train <- xtabs(value~., x_f100_train)
    x_f150_train <- stacked_clr_all_train[stacked_clr_all_train$view == 'f_150',]
    x_f150_train <- subset(x_f150_train, select=c("sample", "feature", "value"))
    x_f150_train <- xtabs(value~., x_f150_train)
    x_f200_train <- stacked_clr_all_train[stacked_clr_all_train$view == 'f_200',]
    x_f200_train <- subset(x_f200_train, select=c("sample", "feature", "value"))
    x_f200_train <- xtabs(value~., x_f200_train)

    # testing
    stacked_clr_all_test <- fread(
      paste('../../data/simulations/lowrank/mofa_tables_fsignal/fold-', fold_,
            '-subset-test-fsignal', nfeats, ".tsv.gz", sep=""), sep="\t")

    # seperate by view
    x_f100_test <- stacked_clr_all_test[stacked_clr_all_test$view == 'f_100',]
    x_f100_test <- subset(x_f100_test, select=c("sample", "feature", "value"))
    x_f100_test <- xtabs(value~., x_f100_test)
    x_f150_test <- stacked_clr_all_test[stacked_clr_all_test$view == 'f_150',]
    x_f150_test <- subset(x_f150_test, select=c("sample", "feature", "value"))
    x_f150_test <- xtabs(value~., x_f150_test)
    x_f200_test <- stacked_clr_all_test[stacked_clr_all_test$view == 'f_200',]
    x_f200_test <- subset(x_f200_test, select=c("sample", "feature", "value"))
    x_f200_test <- xtabs(value~., x_f200_test)

    # compile these into a single X object
    x_f100_test <- as.data.frame.matrix(x_f100_test)
    x_f100_train <- as.data.frame.matrix(x_f100_train)
    x_f100 <- merge(x_f100_test, x_f100_train, all = TRUE)
    x_f100[is.na(x_f100)] <- 0
    rownames(x_f100) <- c(rownames(x_f100_test), rownames(x_f100_train))

    x_f150_test <- as.data.frame.matrix(x_f150_test)
    x_f150_train <- as.data.frame.matrix(x_f150_train)
    x_f150 <- merge(x_f150_test ,x_f150_train, all = TRUE)
    x_f150[is.na(x_f150)] <- 0
    rownames(x_f150) <-  c(rownames(x_f150_test), rownames(x_f150_train))

    x_f200_test <- as.data.frame.matrix(x_f200_test)
    x_f200_train <- as.data.frame.matrix(x_f200_train)
    x_f200 <- merge(x_f200_test, x_f200_train, all = TRUE)
    x_f200[is.na(x_f200)] <- 0
    rownames(x_f200) <-  c(rownames(x_f200_test), rownames(x_f200_train))

    X <- list(f_100 = as.data.frame.matrix(x_f100),
              f_150 = as.data.frame.matrix(x_f150),
              f_200 = as.data.frame.matrix(x_f200))

    # make sure names are shared
    shared_names <- unique(Reduce(intersect,  lapply(X, rownames)))
    X <- list(f_100 = as.data.frame.matrix(x_f100[rownames(x_f100) %in% shared_names, ]),
              f_150 = as.data.frame.matrix(x_f150[rownames(x_f150) %in% shared_names, ]),
              f_200 = as.data.frame.matrix(x_f200[rownames(x_f200) %in% shared_names, ]))

    ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # DIABLO in mixOmics with a value of 0.5 in the design matrix 
    # (cite: https://diabetesjournals.org/care/article/41/10/2178/36703/Intestinal-Metaproteomics-Reveals-Host-Microbiota)
    design = matrix(0.5, 
                    ncol = length(X),
                    nrow = length(X),
                    dimnames = list(names(X),
                                    names(X)))
    diag(design) = 0 # set diagonal to 0s
    model.spls <- block.spls(X, indY = 1, ncomp = 3, design = design)
    factors <- ldply(Map(cbind, model.spls$variates, subject = lapply(model.spls$variates, rownames)))
    write.csv(factors, paste("../../data/simulations/lowrank/mixomics_fsignal/",
                             fold_ ,".factors.model.fsignal", nfeats, ".csv", sep=""))
  }
}
