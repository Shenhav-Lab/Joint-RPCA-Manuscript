# requires numpy==1.23 and scipy==1.8.1 not pinned in the setup.py
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MOFA2")

library(magrittr)
library(reshape2)
library(MOFA2)
library(data.table)

setwd('/Users/bec5786/Desktop/Shenhav Lab/github/Joint-RPCA-Manuscript/code/simulations-benchmarking/')

folds <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10')
fsignals <- c('2', '5', '10', '20', '30', '40', '50', '60')

for (fold_ in folds) {
  print(paste0("Fold: ", fold_))
  
  for (nfeats in fsignals) {
    stacked_clr_all <- fread(
      paste('../../data/simulations/lowrank/mofa_tables_fsignal/fold-',fold_,'-subset-train-fsignal', nfeats,'.tsv.gz', sep=""),
      sep="\t")
    
    # run mofa based on tutorial (https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/microbiome_vignette.html)
    mofa <- create_mofa(stacked_clr_all) 
    print("done")
    model_opts <- get_default_model_options(mofa)
    model_opts$num_factors <- 3
    train_opts <- get_default_training_options(mofa)
    train_opts$convergence_mode <- "medium" #corresponds to the deltaELBO change, medium = 0.00005% - alternatives are "fast" or "slow"
    train_opts$seed <- 42
    mofa <- prepare_mofa(mofa, model_options = model_opts, training_options = train_opts)
    
    t1 <- Sys.time()
    mofa <- run_mofa(mofa, outfile=paste("../../data/simulations/lowrank/mofa_fsignal/",fold_,".model.fsignal", nfeats,".hdf5", sep="",
                                         use_basilisk = TRUE))
    t2 <- Sys.time()
    runtime <- difftime(t2, t1, units = "secs")
    #save runtime
    writeLines(as.character(runtime), paste("../../data/simulations/lowrank/mofa_fsignal/",fold_,".runtime.fsignal", nfeats,".txt", sep=""))
    
    factors <- get_factors(mofa, as.data.frame = TRUE)
    weights <- get_weights(mofa, as.data.frame = TRUE)
    
    write.csv(factors, paste("../../data/simulations/lowrank/mofa_fsignal/",fold_ ,".factors.model.fsignal", nfeats,".csv", sep=""))
    write.csv(weights, paste("../../data/simulations/lowrank/mofa_fsignal/",fold_,".weights.model.fsignal", nfeats,".csv", sep=""))
    }
}
