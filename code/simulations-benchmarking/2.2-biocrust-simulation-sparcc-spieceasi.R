library(phyloseq)
library(igraph)
library(devtools)
library(SpiecEasi)

# dense
metab <- read.delim(
  '../../data/simulations/soil_benchmarking/data/abs.metab.tsv',
  row.names = 1)
microb <- read.delim(
  '../../data/simulations/soil_benchmarking/data/abs.microbes.tsv',
  row.names = 1)

# Fix row names
microb <- microb[,colnames(metab)]

## Load into phyloseq, object used to hold OTU tables
metabolite_phy <- phyloseq(OTU=otu_table(metab, taxa_are_rows=TRUE))
microbe_phy    <- phyloseq(OTU=otu_table(microb, taxa_are_rows=TRUE))

# Sparcc
comb<- rbind(metab,microb)
sparc_test<- sparcc(t(comb),iter=10, th=0)
cor_result<- sparc_test$Cor
rownames(cor_result) <- colnames(cor_result) <- rownames((comb))
write.csv(cor_result, "../../data/simulations/soil_benchmarking/abs/SparCC.csv")

## SPIECEASI
## Load into phyloseq, object used to hold OTU tables
metabolite_phy <- phyloseq(OTU=otu_table(metab, taxa_are_rows=TRUE))
microbe_phy    <- phyloseq(OTU=otu_table(microb, taxa_are_rows=TRUE))
## Run SE on both objects
se.test <- spiec.easi(list(microbe_phy, metabolite_phy), method='mb',
                      nlambda=40, lambda.min.ratio=1e-2,
                      pulsar.params = list(thresh = 0.05, ncores=2))
## Plot using igraph 
dtype <- c(rep(1,ntaxa(microbe_phy)), rep(2,ntaxa(metabolite_phy)))
nodenames <- c(taxa_names(microbe_phy), taxa_names(metabolite_phy))
ig.se <- adj2igraph(getRefit(se.test))
## Get edgelist weights from beta (optimal covariance) matrix 
sebeta <- as.matrix(symBeta(getOptBeta(se.test), mode='maxabs'))
rownames(sebeta) <- colnames(sebeta) <- nodenames
el <- get.edgelist(ig.se)
sizes <- list(rep(1, length(el[,1])))
for (i in 1:length(el[,1])){
  first <- el[,1][i]
  second <- el[,2][i] 
  sizes[i]<-sebeta[first,second]
}
E(ig.se)$weight <- unlist(sizes)
## Write edgelist with weights from beta matrix: model coefficients from neighborhood selection
V(ig.se)$name <- nodenames
ig.el<- as_edgelist(ig.se,names=TRUE)
ig.el.weight <- cbind(ig.el , round(E(ig.se)$weight,5 ))
write.csv(ig.el.weight,"../../data/simulations/soil_benchmarking/abs/SEmultitest.csv")
sebeta <- as.matrix(symBeta(getOptBeta(se.test), mode='maxabs'))
rownames(sebeta) <- colnames(sebeta) <- nodenames
write.csv(sebeta,"../../data/simulations/soil_benchmarking/abs/SPIECEASI.csv")
    

# subsampled
densities <- c('.0.', '.10.')
for (density_ in densities) {
    
    print(density_)
    metab <- read.delim(
        paste("../../data/simulations/soil_benchmarking/data/rel.metabolites", density_,"tsv", sep=""),
      row.names = 1)
    microb <- read.delim(
      paste("../../data/simulations/soil_benchmarking/data/rel.microbes", density_,"tsv", sep=""),
      row.names = 1)
    microb <- subset(microb, rowSums(microb) > 0)
    metab <- subset(metab, rowSums(metab) > 0)
    
    # Fix row names
    microb <- microb[,colnames(metab)]
   
    # Sparcc
    comb<- rbind(metab,microb)
    sparc_test<- sparcc(t(comb), iter=10, th=0)
    cor_result<- sparc_test$Cor
    rownames(cor_result) <- colnames(cor_result) <- rownames((comb))
    write.csv(cor_result, paste("../../data/simulations/soil_benchmarking/rel/SparCC", density_,"csv", sep=""))

    ## SPIECEASI
    ## Load into phyloseq, object used to hold OTU tables
    metabolite_phy <- phyloseq(OTU=otu_table(metab, taxa_are_rows=TRUE))
    microbe_phy    <- phyloseq(OTU=otu_table(microb, taxa_are_rows=TRUE))
    ## Run SE on both objects
    se.test <- spiec.easi(list(microbe_phy, metabolite_phy), 
                          method='mb', nlambda=40, lambda.min.ratio=0.01,
                          pulsar.params = list(thresh = 0.05, ncores=2))
    ## Plot using igraph 
    dtype <- c(rep(1,ntaxa(microbe_phy)), rep(2,ntaxa(metabolite_phy)))
    nodenames <- c(taxa_names(microbe_phy), taxa_names(metabolite_phy))
    ig.se <- adj2igraph(getRefit(se.test))
    ## Get edgelist weights from beta (optimal covariance) matrix 
    sebeta <- as.matrix(symBeta(getOptBeta(se.test), mode='maxabs'))
    rownames(sebeta) <- colnames(sebeta) <- nodenames
    el <- get.edgelist(ig.se)
    sizes <- list(rep(1, length(el[,1])))
    for (i in 1:length(el[,1])){
      first <- el[,1][i]
      second <- el[,2][i] 
      sizes[i]<-sebeta[first,second]
    }
    E(ig.se)$weight <- unlist(sizes)
    ## Write edgelist with weights from beta matrix: model coefficients from neighborhood selection
    V(ig.se)$name <- nodenames
    ig.el<- as_edgelist(ig.se,names=TRUE)
    ig.el.weight <- cbind(ig.el , round(E(ig.se)$weight,5 ))
    write.csv(ig.el.weight, paste("../../data/simulations/soil_benchmarking/rel/SEmultitest", density_,"csv", sep=""))
    sebeta <- as.matrix(symBeta(getOptBeta(se.test), mode='maxabs'))
    rownames(sebeta) <- colnames(sebeta) <- nodenames
    write.csv(sebeta, paste("../../data/simulations/soil_benchmarking/rel/SPIECEASI", density_,"csv", sep=""))

}
