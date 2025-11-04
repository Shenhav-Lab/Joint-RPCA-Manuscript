# requires numpy==1.23 and scipy==1.8.1 not pinned in the setup.py
library(magrittr)
library(reshape2)
library(MOFA2)
library(data.table)

#if running MOFA+ with rclr-transformed table without train/test splits
stacked_combined <- fread(
  '../../data/case-studies/ihmp/mofa-tables/rclr-all-data.tsv.gz', sep="\t")

# run mofa based on tutorial (https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/microbiome_vignette.html)
mofa <- create_mofa(stacked_combined)
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 3
train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode <- "medium"
train_opts$seed <- 42
mofa <- prepare_mofa(mofa, model_options = model_opts, training_options = train_opts)
mofa <- run_mofa(mofa, outfile='../../data/case-studies/ihmp/mofa-results/rclr-alldata-model.hdf5')
factors <- get_factors(mofa, as.data.frame = TRUE)
weights <- get_weights(mofa, as.data.frame = TRUE)
write.csv(factors, '../../data/case-studies/ihmp/mofa-results/rclr-alldata-factors.csv')
write.csv(weights, '../../data/case-studies/ihmp/mofa-results/rclr-alldata-weights.csv')