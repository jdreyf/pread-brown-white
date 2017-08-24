## Author: Jonathan Dreyfuss
## Date updated: 2017-08-24

source("B:/fcns/config.R")

##read
mat <- read.csv("geo_matrix.csv", row.names=1)
annot <- read.csv("probe_annot.csv", row.names=1)
pheno <- read.csv("geo_metadata.csv", row.names=1)

##limma
#paired design
des <- model.matrix(~subject+tissue, data=pheno)
colnames(des) <- gsub("tissue|\\(|\\)", "", colnames(des))
#apply limma
tt <- limma.contrasts(mat, pheno$tissue, contrasts.v = c(WATvsBAT="WAT"), design = des)
tt.df <- data.frame(signif(tt, 3), annot[rownames(tt),])

##Subset to Recon 2.1A genes that are associated with reactions that carry flux in Seahorse media
reduced.g <- read.csv("Recon_21A_reduced_genes.csv")[,1]
reduced.g <- unique(gsub("\\..+", "", reduced.g))

##parse entrez gene id's to see if they match entrez annotation
genes.in.reduced.model <- apply(as.matrix(tt.df$Entrez), 1, FUN=function(x){
  any(unlist(strsplit(x=x, split=" /// ")) %in% reduced.g)
})
#add a column that represents my computation of whether gene is represented in model
tt.df$in.reduced.model <- genes.in.reduced.model
write.csv(tt.df, "pread_gene_stats.csv")

##check
sum(tt.df$WATvsBAT.FDR[tt.df$in.reduced.model] < 0.25)==4
