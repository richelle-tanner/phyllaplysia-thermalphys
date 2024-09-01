#ddRAD analysis: find out the best settings for Stacks
#downloaded genepop files from Comet SDSC XSEDE server (created 9.29.19)
#downloaded M2, M3, M4, M5 - need to figure out best fit before moving forward
#r80 method from: Paris et al. 2017 Methods in Ecology and Evolution
#use these source files: M_populations.sumstats_summary.tsv - seems unnecessary to follow the exact method if all you need is simple plots
#"variant and fixed" table copied into new document: M_r80pop.sumstats_summary.txt

#only changed M, and m was set to match M (per recommendations of Paris et al. 2017)
#We found that the highest polymorphism of r80 loci resulted from n = M, n = M1orn = M + 1(TableS9)andsuggestthisasthebest method for obtaining the optimal value for n.

#looks like M2 is the best, which happens to be the default setting...

#Author: R Tanner
#Original: 11.9.19
#Last Edit: 11.9.19

library(ggplot2)

#import data
ddr <- "~/Documents/R_Data/ddRADSeq/populations_output_r80/"
M2r80 <- read.table(paste(ddr, "M2r80pop.sumstats_summary.txt", sep=''), fill=TRUE, header=TRUE)
M3r80 <- read.table(paste(ddr, "M3r80pop.sumstats_summary.txt", sep=''), fill=TRUE, header=TRUE)
M4r80 <- read.table(paste(ddr, "M4r80pop.sumstats_summary.txt", sep=''), fill=TRUE, header=TRUE)
M5r80 <- read.table(paste(ddr, "M5r80pop.sumstats_summary.txt", sep=''), fill=TRUE, header=TRUE)

#box plot of sites (all: variant and fixed)
sites <- data.frame(pop_ID=M2r80$Pop_ID, setting=rep(c("M2", "M3", "M4", "M5"), each=11),
                      sites=c(M2r80$Sites, M3r80$Sites, M4r80$Sites, M5r80$Sites))

ggplot(sites, aes(x=setting, y=sites)) +
  geom_boxplot() + theme_bw()

#box plot of polymorphic sites (all: variant and fixed)
polym_sites <- data.frame(pop_ID=M2r80$Pop_ID, setting=rep(c("M2", "M3", "M4", "M5"), each=11),
                    polymorphic_sites=c(M2r80$Polymorphic_Sites, M3r80$Polymorphic_Sites, M4r80$Polymorphic_Sites, M5r80$Polymorphic_Sites))

ggplot(polym_sites, aes(x=setting, y=polymorphic_sites)) +
  geom_boxplot() + theme_bw()

#box plot of private sites
private_sites <- data.frame(pop_ID=M2r80$Pop_ID, setting=rep(c("M2", "M3", "M4", "M5"), each=11),
                    private_sites=c(M2r80$Private, M3r80$Private, M4r80$Private, M5r80$Private))

ggplot(private_sites, aes(x=setting, y=private_sites)) +
  geom_boxplot() + theme_bw()
