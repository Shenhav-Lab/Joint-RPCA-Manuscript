library(propr)
library(magrittr)

# dense

met <- read.delim(
  '../../data/simulations/soil_benchmarking/data/abs.metab.tsv',
  row.names = 1)
met <- subset(met, rowSums(met) > 0) %>% t
mic <- read.delim(
  '../../data/simulations/soil_benchmarking/data/abs.microbes.tsv',
  row.names = 1)
mic <- subset(mic, rowSums(mic) > 0) %>% t

# Fix row names
mic <- mic[rownames(met),]

# Replace with smallest non-zero value
met[met == 0] <- min(met[met != 0])
mic[mic == 0] <- min(mic[mic != 0])

# Do a CLR of met and mic separately
clr <- function(x) sweep(log(x), 1, rowMeans(log(x)), "-")
agg <- cbind(clr(met), clr(mic))

# PHI
pr <- propr:::lr2phi(agg)
colnames(pr) <- colnames(agg)
rownames(pr) <- colnames(agg)
write.csv(pr, "../../data/simulations/soil_benchmarking/abs/prop_matrix_PHI.csv")

# RHO
pr <- propr:::lr2rho(agg)
colnames(pr) <- colnames(agg)
rownames(pr) <- colnames(agg)
write.csv(pr, "../../data/simulations/soil_benchmarking/abs/prop_matrix_RHO.csv")


# PEARSON
pr <- stats::cor(agg)
colnames(pr) <- colnames(agg)
write.csv(pr, "../../data/simulations/soil_benchmarking/abs/prop_matrix_pearson.csv")

# SPEARMAN
pr <- stats::cor(agg, method = "spearman")
colnames(pr) <- colnames(agg)
write.csv(pr, "../../data/simulations/soil_benchmarking/abs/prop_matrix_spearman.csv")

# subsampled
densities <- c('.0.', '.10.')
for (density_ in densities) {
    
    print(density_)
    met <- read.delim(
        paste("../../data/simulations/soil_benchmarking/data/rel.metabolites", density_,"tsv", sep=""),
      row.names = 1)
    met <- subset(met, rowSums(met) > 0) %>% t
    mic <- read.delim(
      paste("../../data/simulations/soil_benchmarking/data/rel.microbes", density_,"tsv", sep=""),
      row.names = 1)
    mic <- subset(mic, rowSums(mic) > 0) %>% t

    # Fix row names
    mic <- mic[rownames(met),]

    # Replace with smallest non-zero value
    met[met == 0] <- min(met[met != 0])
    mic[mic == 0] <- min(mic[mic != 0])

    # Do a CLR of met and mic separately
    clr <- function(x) sweep(log(x), 1, rowMeans(log(x)), "-")
    agg <- cbind(clr(met), clr(mic))

    # PHI
    pr <- propr:::lr2phi(agg)
    colnames(pr) <- colnames(agg)
    rownames(pr) <- colnames(agg)
    write.csv(pr, paste("../../data/simulations/soil_benchmarking/rel/prop_matrix_PHI", density_,"csv", sep=""))

    # RHO
    pr <- propr:::lr2rho(agg)
    colnames(pr) <- colnames(agg)
    rownames(pr) <- colnames(agg)
    write.csv(pr, paste("../../data/simulations/soil_benchmarking/rel/prop_matrix_RHO", density_,"csv", sep=""))


    # PEARSON
    pr <- stats::cor(agg)
    colnames(pr) <- colnames(agg)
    write.csv(pr, paste("../../data/simulations/soil_benchmarking/rel/prop_matrix_pearson", density_,"csv", sep=""))

    # SPEARMAN
    pr <- stats::cor(agg, method = "spearman")
    colnames(pr) <- colnames(agg)
    write.csv(pr, paste("../../data/simulations/soil_benchmarking/rel/prop_matrix_spearman", density_,"csv", sep=""))

}
