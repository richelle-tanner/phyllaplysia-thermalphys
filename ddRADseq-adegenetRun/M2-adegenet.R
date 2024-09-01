#ddRAD analysis: adegenet
#downloaded genepop files from Comet SDSC XSEDE server (created 11.7.19)
#settings from Stacks "populations" script: 
#final decision: using M2

#from Stacks page, how to load:
#first rename genepop file with gen suffix
#library("adegenet")
#x <- read.genepop("populations.gen")

#Author: R Tanner 
#Original: 11.9.19
#Last Edit: 11.9.19

#load packages
.libPaths("/home/richelle.tanner/scripts/ddRAD/lib/R_libs/")
library(adegenet)

#load data
ddr <- "/home/richelle.tanner/scripts/ddRAD/"
snps <- read.genepop(paste(ddr,"M2populations.snps.gen", sep=''))

#follow the tutorial for PCA: pg 51 of tutorial-basics.pdf "summary of the genetic diversity among the sampled individuals"
#replace missing data by the mean allele frequency
X <- tab(snps, freq=TRUE, NA.method = "mean")
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)

pdf('barplot_pca.pdf')
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))
s.label(pca1$li)
title("PCA of Latitudinal dataset\naxes 1-2") 
add.scatter.eig(pca1$eig[1:20], 3,1,2)
dev.off()

pdf('full_pca.pdf')
s.class(pca1$li, pop(snps))
title("PCA of Latitudinal dataset\naxes 1-2") 
add.scatter.eig(pca1$eig[1:20], 3,1,2)
dev.off()

table <- snps@tab
write.table(table, "allele_freq_table.txt", quote=F, row.names=F, sep='\t')
