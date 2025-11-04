# requires numpy==1.23 and scipy==1.8.1 not pinned in the setup.py
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MOFA2")

library(magrittr)
library(reshape2)
library(MOFA2)
library(data.table)

setwd('/Users/bec5786/Desktop/Shenhav Lab/Joint-RPCA/joint-rpca-bencmarks/code/simulations-benchmarking/')

folds <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10')
for (fold_ in folds) {
  print(paste0("Fold: ", fold_))
  densities <- c('','3','5','7','9','11')
for (density_ in densities) {
  # dense
  stacked_clr_all <- fread(
    paste('../../data/simulations/ihmp/mofa-reformatted/fold-', fold_,'-subset-train-meta_g_taxonomic_profiles_', density_, "-test.tsv.gz", sep=""),
    sep="\t")
  # run mofa based on tutorial (https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/microbiome_vignette.html)
  mofa <- create_mofa(stacked_clr_all)
  print("done")
  model_opts <- get_default_model_options(mofa)
  model_opts$num_factors <- 3
  train_opts <- get_default_training_options(mofa)
  train_opts$convergence_mode <- "medium"
  train_opts$seed <- 42
  mofa <- prepare_mofa(mofa, model_options = model_opts, training_options = train_opts)
  mofa <- run_mofa(mofa, outfile=paste("../../results/simulations/ihmp/mofa/", fold_, ".model.", density_, ".hdf5", sep="",
                                       use_basilisk = TRUE))
  factors <- get_factors(mofa, as.data.frame = TRUE)
  weights <- get_weights(mofa, as.data.frame = TRUE)
  write.csv(factors, paste("../../results/simulations/ihmp/mofa/", fold_ ,".factors.model.", density_, ".csv", sep=""))
  write.csv(weights, paste("../../results/simulations/ihmp/mofa/", fold_,".weights.model.", density_, ".csv", sep=""))

  }
}
