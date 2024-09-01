#try the 70% file from the Stacks Pipeline (2nd one, for filtering)

library("adegenet")
library("stringr")
library("ade4")
library("poppr")
library("ape") # To visualize the tree using the "nj" function
library("magrittr")
library("LDcorSV")
library("pegas")
library("hierfstat")

####import data####
ddr <- "~/Documents/R_Data/ddRADSeq/"
#use structure file exported by filtering pipeline
snps <- read.structure(paste(ddr, "populations_r100_m70_firstSNP.str", sep=''))
#complete call was : n.ind = 109, n.loc = 80344, col.pop = 0, col.lab = 1, row.marknames = 0, onerowperind = FALSE
#used Filter_All_tsv.summary.m70_firstSNP.log for individual and loci info
pop_map <- read.delim(paste(ddr, "popmap_16_17/full_popmap_2016_2017.txt", sep=''), header=FALSE)
pop_map_SF <- read.delim(paste(ddr, "popmap_16_17/full_popmapSF_2016_2017.txt", sep=''), header=FALSE)
phenotypes <- read.delim(paste(ddr, "combined_respRateCTmax.txt", sep=''), header=TRUE, sep='\t') #there are spaces in the location names, so always use tab separator
matchID <- read.delim(paste(ddr, "matchIDs_genotype_phenotype.txt", sep=''), header=TRUE, sep='\t') #there are a ton of empty rows at the end

####data cleaning####
#assign populations to the "snps" genind object
included_ind <- as.character(snps$pop)
pop_map_red <- pop_map[which(pop_map$V1 %in% included_ind),]
pop_mapSF_red <- pop_map_SF[which(pop_map_SF$V1 %in% included_ind),]
#also use strata so we can take out the 2016 ones only
pop_map_red$V2 <- as.character(pop_map_red$V2)
pop_mapSF_red$V2 <- as.character(pop_mapSF_red$V2)
strata_df <- data.frame(Population=substr(pop_map_red$V2,1,nchar(pop_map_red$V2)-4), 
                        Year=str_sub(pop_map_red$V2, start=-4),
                        PopSF=substr(pop_mapSF_red$V2,1,nchar(pop_mapSF_red$V2)-4))
strata(snps) <- strata_df #this is a data frame with columns for population and year
nameStrata(snps) <- ~Population/Year/PopSF

match_snpsPOP <- data.frame(ID=included_ind, strata_df)
match2016_snpsPOP <- subset(match_snpsPOP, Year == "2016")

#take out only 2016 - fewer than individuals with phenotypic data (88 here vs 173 phenotypes)
snps_sep <- seppop(snps, ~Year) #subset using strata
snps2016 <- snps_sep$`2016`
#Jennifer's suggestion to scale and clean data before generating a matrix for PCA plot
snps2016_scale <- scaleGen(snps2016, NA.method = "mean") #removes or fills in (with mean) missing alleles, this outputs something that looks like a @tab object
dim(snps2016_scale)
#88 124400

colnames(matchID) <- c("genID", "treatment") #to match the labels in phenotype data frame
matchID$treatment <- as.character(matchID$treatment)
phenotypes$treatment <- as.character(phenotypes$treatment)
phenotypes_genID <- merge(phenotypes, matchID, by="treatment") #there are 173 left - this was cross checked against the library organization list 3.19.20

#match snps2016_scale to the 173 individuals we have phenotypic data for
#this needs to be done before the individual names are lost?
phenotypes_IDs_allAcclim <- as.character(phenotypes_genID$genID)
#get rid of spaces so they match
searchString <- ' '
replacementString <- ''
phenotypes_IDs_allAcclim <- sub(searchString, replacementString, phenotypes_IDs_allAcclim)
genind_names <- as.character(snps@pop)
select <- which(genind_names %in% phenotypes_IDs_allAcclim) #why are there only 72 here? this is because we lost individuals to filtering!
match_snps <- snps[select,]
match_snps_scale <- scaleGen(match_snps, NA.method = "mean") #108830 loci
#also have to reduce phenotype data to available genotypes
phenotypes_genID16 <- phenotypes_genID[which(phenotypes_IDs_allAcclim %in% match2016_snpsPOP$ID),] #this works
#reassign so downstream still works
phenotypes_genID <- phenotypes_genID16

#also find indices separated by acclimation
phenotypes_genID_13 <- subset(phenotypes_genID, temp == 13) #23
phenotypes_genID_17 <- subset(phenotypes_genID, temp == 17) #27
phenotypes_genID_21 <- subset(phenotypes_genID, temp == 21) #22
phenotypes_IDs_13 <- as.character(phenotypes_genID_13$genID)
phenotypes_IDs_17 <- as.character(phenotypes_genID_17$genID)
phenotypes_IDs_21 <- as.character(phenotypes_genID_21$genID)
phenotypes_IDs_13 <- sub(searchString, replacementString, phenotypes_IDs_13)
phenotypes_IDs_17 <- sub(searchString, replacementString, phenotypes_IDs_17)
phenotypes_IDs_21 <- sub(searchString, replacementString, phenotypes_IDs_21)
select13 <- which(genind_names %in% phenotypes_IDs_13) #42
select17 <- which(genind_names %in% phenotypes_IDs_17) #51
select21 <- which(genind_names %in% phenotypes_IDs_21) #43
match_snps13 <- snps[select13,]
match_snps13_scale <- scaleGen(match_snps13, NA.method = "mean") #42724 loci
match_snps17 <- snps[select17,]
match_snps17_scale <- scaleGen(match_snps17, NA.method = "mean") #55786 loci
match_snps21 <- snps[select21,]
match_snps21_scale <- scaleGen(match_snps21, NA.method = "mean") #48826 loci

#cut down to remaining genotypes for phenotypic data
match_snps13_names <- genind_names[which(genind_names %in% phenotypes_IDs_13)]
match_snps17_names <- genind_names[which(genind_names %in% phenotypes_IDs_17)]
match_snps21_names <- genind_names[which(genind_names %in% phenotypes_IDs_21)]
phenotypes_genID_13$genID <- sub(searchString, replacementString, phenotypes_genID_13$genID)
phenotypes_genID_17$genID <- sub(searchString, replacementString, phenotypes_genID_17$genID)
phenotypes_genID_21$genID <- sub(searchString, replacementString, phenotypes_genID_21$genID)
phenotypes_genID_select13 <- phenotypes_genID_13[which(phenotypes_genID_13$genID %in% match_snps13_names),]
phenotypes_genID_select17 <- phenotypes_genID_17[which(phenotypes_genID_17$genID %in% match_snps17_names),]
phenotypes_genID_select21 <- phenotypes_genID_21[which(phenotypes_genID_21$genID %in% match_snps21_names),]
phenotypes_genID_select13 <- phenotypes_genID_select13[order(phenotypes_genID_select13$genID),]
phenotypes_genID_select17 <- phenotypes_genID_select17[order(phenotypes_genID_select17$genID),]
phenotypes_genID_select21 <- phenotypes_genID_select21[order(phenotypes_genID_select21$genID),]


match_snps_names <- genind_names[which(genind_names %in% phenotypes_IDs_allAcclim)]
phenotypes_genID$genID <- sub(searchString, replacementString, phenotypes_genID$genID)
phenotypes_genID_select <- phenotypes_genID[which(phenotypes_genID$genID %in% match_snps_names),]
#make sure the phenotypic data is in the same order as the genotypic data
phenotypes_genID_select <- phenotypes_genID_select[order(phenotypes_genID_select$genID),]

####generate distance matrices for each dataset####

#using the output from Jennifer's data type that is already scaled
distPop <- dist(match_snps_scale, method="euclidean")

#generate PCA for Rauri
genind_names_phen <- genind_names[which(genind_names %in% phenotypes_IDs_allAcclim)]
genPop_names <- t(data.frame(strsplit(genind_names_phen, "_")))
genPop_names <- unname(genPop_names[,1])
mod <- betadisper(distPop, genPop_names, type="median")
plot(mod)
boxplot(mod)

#phenotype
#scale data first!
#we want to keep respRate_noHS, change, CTmax, weight_mg (this is the order of the columns since the names go away with the matrix conversion)
#do we need to correct somehow for acclimation temperature? or at least do separate comparisons
phenotypes_allAcclim <- as.matrix(phenotypes_genID_select[,c(4,6,9,11)])
phenotypes_allAcclim <- scale(phenotypes_allAcclim)
distTrait <- dist(phenotypes_allAcclim, method="euclidean")

#what if you do the phenotypes separately
respRate_allAcclim <- phenotypes_genID_select[,4]
respRate_allAcclim <- scale(respRate_allAcclim)
distRespRate <- dist(respRate_allAcclim, method="euclidean")
respChange_allAcclim <- phenotypes_genID_select[,6]
respChange_allAcclim <- scale(respChange_allAcclim)
distRespChange <- dist(respChange_allAcclim, method="euclidean")
CTmax_allAcclim <- phenotypes_genID_select[,9]
CTmax_allAcclim <- scale(CTmax_allAcclim)
distCTmax <- dist(CTmax_allAcclim, method="euclidean")

####treatment-specific distance matrices####
distPop13 <- dist(match_snps13_scale, method="euclidean")
distPop17 <- dist(match_snps17_scale, method="euclidean")
distPop21 <- dist(match_snps21_scale, method="euclidean")

phenotypes_13 <- as.matrix(phenotypes_genID_select13[,c(4,6,9,11)])
phenotypes_13 <- scale(phenotypes_13)
distTrait13 <- dist(phenotypes_13, method="euclidean")

#what if you do the phenotypes separately
respRate_13 <- phenotypes_genID_select13[,4]
respRate_13 <- scale(respRate_13)
distRespRate13 <- dist(respRate_13, method="euclidean")
respChange_13 <- phenotypes_genID_select13[,6]
respChange_13 <- scale(respChange_13)
distRespChange13 <- dist(respChange_13, method="euclidean")
CTmax_13 <- phenotypes_genID_select13[,9]
CTmax_13 <- scale(CTmax_13)
distCTmax13 <- dist(CTmax_13, method="euclidean")

phenotypes_17 <- as.matrix(phenotypes_genID_select17[,c(4,6,9,11)])
phenotypes_17 <- scale(phenotypes_17)
distTrait17 <- dist(phenotypes_17, method="euclidean")

#what if you do the phenotypes separately
respRate_17 <- phenotypes_genID_select17[,4]
respRate_17 <- scale(respRate_17)
distRespRate17 <- dist(respRate_17, method="euclidean")
respChange_17 <- phenotypes_genID_select17[,6]
respChange_17 <- scale(respChange_17)
distRespChange17 <- dist(respChange_17, method="euclidean")
CTmax_17 <- phenotypes_genID_select17[,9]
CTmax_17 <- scale(CTmax_17)
distCTmax17 <- dist(CTmax_17, method="euclidean")

phenotypes_21 <- as.matrix(phenotypes_genID_select21[,c(4,6,9,11)])
phenotypes_21 <- scale(phenotypes_21)
distTrait21 <- dist(phenotypes_21, method="euclidean")

#what if you do the phenotypes separately
respRate_21 <- phenotypes_genID_select21[,4]
respRate_21 <- scale(respRate_21)
distRespRate21 <- dist(respRate_21, method="euclidean")
respChange_21 <- phenotypes_genID_select21[,6]
respChange_21 <- scale(respChange_21)
distRespChange21 <- dist(respChange_21, method="euclidean")
CTmax_21 <- phenotypes_genID_select21[,9]
CTmax_21 <- scale(CTmax_21)
distCTmax21 <- dist(CTmax_21, method="euclidean")


####mantel.rtest####
#look at ddRADgen_phen_results - Mantel test results
#make sure the distance matrices are in the same order for individuals
#cursory look shows that there is no correlation between any of the traits (maybe respiration change?) and genetics
#perhaps we need to scale these values based on acclimation condition to see anything?
#what about summarizing by population?
r1 <- mantel.rtest(distPop, distRespRate, nrepet = 999)
plot(r1, main = "Mantel's test comparing\ngenetics and multivariate\nthermal tolerance phenotype")

#treatment-specific (separated by acclimation)
#still nothing significant!
r2 <- mantel.rtest(distPop13, distRespChange13, nrepet = 999)
plot(r2, main = "Mantel's test comparing\ngenetics and HS RR trait thermal\ntolerance phenotype in 13C acclim")

r3 <- mantel.rtest(distPop13, distRespRate13, nrepet = 999)
plot(r3, main = "Mantel's test comparing\ngenetics and RR trait thermal\ntolerance phenotype in 13C acclim")

r4 <- mantel.rtest(distPop13, distCTmax13, nrepet = 999)
plot(r4, main = "Mantel's test comparing\ngenetics and CTmax trait thermal\ntolerance phenotype in 13C acclim")

r5 <- mantel.rtest(distPop17, distRespChange17, nrepet = 999)
plot(r5, main = "Mantel's test comparing\ngenetics and HS RR trait thermal\ntolerance phenotype in 17C acclim")

r6 <- mantel.rtest(distPop17, distRespRate17, nrepet = 999)
plot(r6, main = "Mantel's test comparing\ngenetics and RR trait thermal\ntolerance phenotype in 17C acclim")

r7 <- mantel.rtest(distPop17, distCTmax17, nrepet = 999)
plot(r7, main = "Mantel's test comparing\ngenetics and CTmax trait thermal\ntolerance phenotype in 17C acclim")

r8 <- mantel.rtest(distPop21, distRespChange21, nrepet = 999) #SIGNIFICANT
plot(r8, main = "Mantel's test comparing\ngenetics and HS RR trait thermal\ntolerance phenotype in 21C acclim")

r9 <- mantel.rtest(distPop21, distRespRate21, nrepet = 999)
plot(r9, main = "Mantel's test comparing\ngenetics and RR trait thermal\ntolerance phenotype in 21C acclim")

r10 <- mantel.rtest(distPop21, distCTmax21, nrepet = 999)
plot(r10, main = "Mantel's test comparing\ngenetics and CTmax trait thermal\ntolerance phenotype in 21C acclim")

####inbreeding####
#separate out populations
snps_sep2016 <- seppop(snps2016, ~Population) #subset using strata
snps_sep2016SF <- seppop(snps2016, ~PopSF) #so we can combine the SF population
SF16 <- snps_sep2016SF$`SF`
Albany16 <- snps_sep2016$`Albany`
BayFarm16 <- snps_sep2016$`BayFarm`
Coos16 <- snps_sep2016$`Coos`
CrownBeach16 <- snps_sep2016$`CrownBeach`
Humboldt16 <- snps_sep2016$`Humboldt`
Morro16 <- snps_sep2016$`Morro`
OceanShores16 <- snps_sep2016$`OceanShores`
Paradise16 <- snps_sep2016$`Paradise`
PtMolate16 <- snps_sep2016$`PointMolate`
PtOrient16 <- snps_sep2016$`PointOrient`
Tomales16 <- snps_sep2016$`Tomales`

#SF 16 combined
SF16_in <- inbreeding(SF16, N=100)
head(SF16_in[[1]],20)
Fbar_SF16 <- sapply(SF16_in, mean)
Fbar_SF16 <- data.frame(pop=rep("SF", length(Fbar_SF16)), inbreed=Fbar_SF16)
hist(Fbar_SF16$inbreed, col="firebrick", main="Average inbreeding in \n2016 SF population")
mean(Fbar_SF16$inbreed) #0.5027754
#SF16_F <- Fst(as.loci(SF16)) #SAYS THEY ARE NOT DIPLOID BUT THEY ARE??

#Albany 16
Alb16_in <- inbreeding(Albany16, N=100)
head(Alb16_in[[1]],20)
Fbar_Alb16 <- sapply(Alb16_in, mean)
Fbar_Alb16 <- data.frame(pop=rep("Albany", length(Fbar_Alb16)), inbreed=Fbar_Alb16)
hist(Fbar_Alb16$inbreed, col="firebrick", main="Average inbreeding in \n2016 Albany population")
mean(Fbar_Alb16$inbreed) #0.5111653
#Alb16_F <- Fst(as.loci(Albany16redL))

#BayFarm 16
BF16_in <- inbreeding(BayFarm16, N=100)
head(BF16_in[[1]],20)
Fbar_BF16 <- sapply(BF16_in, mean)
Fbar_BF16 <- data.frame(pop=rep("BayFarm", length(Fbar_BF16)), inbreed=Fbar_BF16)
hist(Fbar_BF16$inbreed, col="firebrick", main="Average inbreeding in \n2016 BayFarm population")
mean(Fbar_BF16$inbreed) #0.4989537
#BF16_F <- Fst(as.loci(BayFarm16redL))

#Coos 16
Coo16_in <- inbreeding(Coos16, N=100)
head(Coo16_in[[1]],20)
Fbar_Coo16 <- sapply(Coo16_in, mean)
Fbar_Coo16 <- data.frame(pop=rep("Coos", length(Fbar_Coo16)), inbreed=Fbar_Coo16)
hist(Fbar_Coo16$inbreed, col="firebrick", main="Average inbreeding in \n2016 Coos population")
mean(Fbar_Coo16$inbreed) #0.4938794
#Coo16_F <- Fst(as.loci(Coos16))

#CrownBeach 16
CB16_in <- inbreeding(CrownBeach16, N=100)
head(CB16_in[[1]],20)
Fbar_CB16 <- sapply(CB16_in, mean)
Fbar_CB16 <- data.frame(pop=rep("CrownBeach", length(Fbar_CB16)), inbreed=Fbar_CB16)
hist(Fbar_CB16$inbreed, col="firebrick", main="Average inbreeding in \n2016 CrownBeach population")
mean(Fbar_CB16$inbreed) #0.5043298
#CB16_F <- Fst(as.loci(CrownBeach16))

#Humboldt 16
Hum16_in <- inbreeding(Humboldt16, N=100)
head(Hum16_in[[1]],20)
Fbar_Hum16 <- sapply(Hum16_in, mean)
Fbar_Hum16 <- data.frame(pop=rep("Humboldt", length(Fbar_Hum16)), inbreed=Fbar_Hum16)
hist(Fbar_Hum16$inbreed, col="firebrick", main="Average inbreeding in \n2016 Humboldt population")
mean(Fbar_Hum16$inbreed) #0.5235243
#Hum16_F <- Fst(as.loci(Humboldt16))

#Morro 16
Mor16_in <- inbreeding(Morro16, N=100)
head(Mor16_in[[1]],20)
Fbar_Mor16 <- sapply(Mor16_in, mean)
Fbar_Mor16 <- data.frame(pop=rep("Morro", length(Fbar_Mor16)), inbreed=Fbar_Mor16)
hist(Fbar_Mor16$inbreed, col="firebrick", main="Average inbreeding in \n2016 Morro population")
mean(Fbar_Mor16$inbreed) #0.5202453
#Mor16_F <- Fst(as.loci(Morro16))

#OceanShores 16
OS16_in <- inbreeding(OceanShores16, N=100)
head(OS16_in[[1]],20)
Fbar_OS16 <- sapply(OS16_in, mean)
Fbar_OS16 <- data.frame(pop=rep("OceanShores", length(Fbar_OS16)), inbreed=Fbar_OS16)
hist(Fbar_OS16$inbreed, col="firebrick", main="Average inbreeding in \n2016 OceanShores population")
mean(Fbar_OS16$inbreed) #0.5218347
#OS16_F <- Fst(as.loci(OceanShores16))

#Paradise 16
Par16_in <- inbreeding(Paradise16, N=100)
head(Par16_in[[1]],20)
Fbar_Par16 <- sapply(Par16_in, mean)
Fbar_Par16 <- data.frame(pop=rep("Paradise", length(Fbar_Par16)), inbreed=Fbar_Par16)
hist(Fbar_Par16$inbreed, col="firebrick", main="Average inbreeding in \n2016 Paradise population")
mean(Fbar_Par16$inbreed) #0.495499
#Par16_F <- Fst(as.loci(Paradise16))

#PtMolate 16
PM16_in <- inbreeding(PtMolate16, N=100)
head(PM16_in[[1]],20)
Fbar_PM16 <- sapply(PM16_in, mean)
Fbar_PM16 <- data.frame(pop=rep("PtMolate", length(Fbar_PM16)), inbreed=Fbar_PM16)
hist(Fbar_PM16$inbreed, col="firebrick", main="Average inbreeding in \n2016 PtMolate population")
mean(Fbar_PM16$inbreed) #0.5002837
#PM16_F <- Fst(as.loci(PointMolate16))

#PtOrient 16
PO16_in <- inbreeding(PtOrient16, N=100)
head(PO16_in[[1]],20)
Fbar_PO16 <- sapply(PO16_in, mean)
Fbar_PO16 <- data.frame(pop=rep("PtOrient", length(Fbar_PO16)), inbreed=Fbar_PO16)
hist(Fbar_PO16$inbreed, col="firebrick", main="Average inbreeding in \n2016 PtOrient population")
mean(Fbar_PO16$inbreed) #0.495103
#PO16_F <- Fst(as.loci(PtOrient16))

#Tomales 16
Tom16_in <- inbreeding(Tomales16, N=100)
head(Tom16_in[[1]],20)
Fbar_Tom16 <- sapply(Tom16_in, mean)
Fbar_Tom16 <- data.frame(pop=rep("Tomales", length(Fbar_Tom16)), inbreed=Fbar_Tom16)
hist(Fbar_Tom16$inbreed, col="firebrick", main="Average inbreeding in \n2016 Tomales population")
mean(Fbar_Tom16$inbreed) #0.479597
#Tom16_F <- Fst(as.loci(Tomales16))

##compare averages and variances for inbreeding between populations##
#make DF of Fbar results (summarizes across loci for each individual)
all_inbreed <- rbind(Fbar_Alb16, Fbar_BF16, Fbar_Coo16, Fbar_CB16, Fbar_Hum16, Fbar_Mor16, Fbar_OS16, Fbar_Par16, Fbar_PM16, Fbar_PO16, Fbar_Tom16)
all_inbreedSF <- rbind(Fbar_SF16, Fbar_Coo16, Fbar_Hum16, Fbar_Mor16, Fbar_OS16, Fbar_Tom16)

#calculate residuals so we can do an anova on the variance as well
Fbar_Alb16_med <- data.frame(pop=Fbar_Alb16$pop, inbreed_residual=abs(Fbar_Alb16$inbreed - median(Fbar_Alb16$inbreed)))
Fbar_BF16_med <- data.frame(pop=Fbar_BF16$pop, inbreed_residual=abs(Fbar_BF16$inbreed - median(Fbar_BF16$inbreed)))
Fbar_Coo16_med <- data.frame(pop=Fbar_Coo16$pop, inbreed_residual=abs(Fbar_Coo16$inbreed - median(Fbar_Coo16$inbreed)))
Fbar_CB16_med <- data.frame(pop=Fbar_CB16$pop, inbreed_residual=abs(Fbar_CB16$inbreed - median(Fbar_CB16$inbreed)))
Fbar_Hum16_med <- data.frame(pop=Fbar_Hum16$pop, inbreed_residual=abs(Fbar_Hum16$inbreed - median(Fbar_Hum16$inbreed)))
Fbar_Mor16_med <- data.frame(pop=Fbar_Mor16$pop, inbreed_residual=abs(Fbar_Mor16$inbreed - median(Fbar_Mor16$inbreed)))
Fbar_OS16_med <- data.frame(pop=Fbar_OS16$pop, inbreed_residual=abs(Fbar_OS16$inbreed - median(Fbar_OS16$inbreed)))
Fbar_Par16_med <- data.frame(pop=Fbar_Par16$pop, inbreed_residual=abs(Fbar_Par16$inbreed - median(Fbar_Par16$inbreed)))
Fbar_PM16_med <- data.frame(pop=Fbar_PM16$pop, inbreed_residual=abs(Fbar_PM16$inbreed - median(Fbar_PM16$inbreed)))
Fbar_PO16_med <- data.frame(pop=Fbar_PO16$pop, inbreed_residual=abs(Fbar_PO16$inbreed - median(Fbar_PO16$inbreed)))
Fbar_Tom16_med <- data.frame(pop=Fbar_Tom16$pop, inbreed_residual=abs(Fbar_Tom16$inbreed - median(Fbar_Tom16$inbreed)))
Fbar_SF16_med <- data.frame(pop=Fbar_SF16$pop, inbreed_residual=abs(Fbar_SF16$inbreed - median(Fbar_SF16$inbreed)))

all_inbreed_res <- rbind(Fbar_Alb16_med, Fbar_BF16_med, Fbar_Coo16_med, Fbar_CB16_med, Fbar_Hum16_med, Fbar_Mor16_med, Fbar_OS16_med, Fbar_Par16_med, Fbar_PM16_med, Fbar_PO16_med, Fbar_Tom16_med)
all_inbreedSF_res <- rbind(Fbar_SF16_med, Fbar_Coo16_med, Fbar_Hum16_med, Fbar_Mor16_med, Fbar_OS16_med, Fbar_Tom16_med)

#no datasets are significant
summary(aov(inbreed~pop, all_inbreed))
summary(aov(inbreed_residual~pop, all_inbreed_res))
summary(aov(inbreed~pop, all_inbreedSF))
summary(aov(inbreed_residual~pop, all_inbreedSF_res))

#export inbreeding coefficients and residuals
inbreed_export <- data.frame(pop=all_inbreed$pop, popSF=all_inbreedSF$pop, 
                             inbreed=all_inbreed$inbreed, inbreedSF=all_inbreedSF$inbreed,
                             inbreed_res=all_inbreed_res$inbreed_residual, inbreedSF_res=all_inbreedSF_res$inbreed_residual)
write.table(inbreed_export, "~/Documents/R_Data/ddRADSeq/inbreeding_coefficients_m70_6.11.20.txt", sep='\t', quote=F, row.names = F)

####private alleles####
#tutorial: https://rdrr.io/cran/poppr/man/private_alleles.html 
private16_default <- private_alleles(snps2016, locus ~ Population, count.alleles = TRUE) #this is supposed to be observed alleles per locus
private16 <- private_alleles(snps2016, locus ~ Population, count.alleles = FALSE) #this was recommended- but why? says this is the raw number of private alleles per locus
#sum number of private alleles per population
total_private16 <- rowSums(private16)
#get percentages - not sure why we want this?
private16sum <- sweep(private16, 2, nAll(snps2016)[colnames(private16)], FUN = "/")

#also do it with SF pooled - use PopSF
private16SF <- private_alleles(snps2016, locus ~ PopSF, count.alleles = FALSE)
total_private16SF <- rowSums(private16SF)

#are there significantly more private alleles in SF groups than outside SF?
private_df <- data.frame(pop=names(total_private16), 
                         SF=c("SF", "SF", "noSF", "SF", "noSF", "noSF", "noSF", "SF", "SF", "SF", "noSF"),
                         alleles=unname(total_private16))
summary(aov(alleles ~ SF, private_df))

####remove missing data####
#https://grunwaldlab.github.io/poppr/reference/missingno.html
#this also might help with errors later on
snps2016redLoci <- missingno(snps2016, "loci") #removes all loci containing missing data, this cuts down to 88 loci/10920 alleles
snps2016redMean <- missingno(snps2016, "mean") #replaces all NAs with mean of alleles for entire dataset - Replaced 7409176 missing values.
#other options are "genotype" (removes individuals with missing data, don't want this - we already know how much the dataset is reduced when you do this) and "zero" (don't want more artificial diversity)

#separate out populations
snps_sep2016redMean <- seppop(snps2016redMean, ~Population) #subset using strata
snps_sep2016SFredMean <- seppop(snps2016redMean, ~PopSF) #so we can combine the SF population
SF16red <- snps_sep2016SFredMean$`SF`
Albany16red <- snps_sep2016redMean$`Albany`
BayFarm16red <- snps_sep2016redMean$`BayFarm`
Coos16red <- snps_sep2016redMean$`Coos`
CrownBeach16red <- snps_sep2016redMean$`CrownBeach`
Humboldt16red <- snps_sep2016redMean$`Humboldt`
Morro16red <- snps_sep2016redMean$`Morro`
OceanShores16red <- snps_sep2016redMean$`OceanShores`
Paradise16red <- snps_sep2016redMean$`Paradise`
PtMolate16red <- snps_sep2016redMean$`PointMolate`
PtOrient16red <- snps_sep2016redMean$`PointOrient`
Tomales16red <- snps_sep2016redMean$`Tomales`

snps_sep2016redLoci <- seppop(snps2016redLoci, ~Population) #subset using strata
snps_sep2016SFredLoci <- seppop(snps2016redLoci, ~PopSF) #so we can combine the SF population
SF16redL <- snps_sep2016SFredLoci$`SF`
Albany16redL <- snps_sep2016redLoci$`Albany`
BayFarm16redL <- snps_sep2016redLoci$`BayFarm`
Coos16redL <- snps_sep2016redLoci$`Coos`
CrownBeach16redL <- snps_sep2016redLoci$`CrownBeach`
Humboldt16redL <- snps_sep2016redLoci$`Humboldt`
Morro16redL <- snps_sep2016redLoci$`Morro`
OceanShores16redL <- snps_sep2016redLoci$`OceanShores`
Paradise16redL <- snps_sep2016redLoci$`Paradise`
PtMolate16redL <- snps_sep2016redLoci$`PointMolate`
PtOrient16redL <- snps_sep2016redLoci$`PointOrient`
Tomales16redL <- snps_sep2016redLoci$`Tomales`

####AMOVA####
#uses poppr package
#https://grunwaldlab.github.io/poppr/reference/poppr.amova.html 
amova.result <- poppr.amova(snps2016, ~Population)
#Found 16870928 missing values.
#53109 loci contained missing values greater than 5%
amova.test <- randtest(amova.result, nrepet = 999)

amova.result.redLoci <- poppr.amova(snps2016redLoci, ~Population)
amova.result.redMean <- poppr.amova(snps2016redMean, ~Population, within=FALSE) #the error says the estimates are skewed and to use the "within=FALSE" option

amova.test.redLoci <- randtest(amova.result.redLoci)
amova.test.redMean <- randtest(amova.result.redMean)

####obs. het.####
#http://adegenet.r-forge.r-project.org/files/montpellier/practical-MVAintro.1.0.pdf

snps2016@pop <- snps2016@strata$Population

full_hwTest <- gstat.randtest(snps2016)

#this uses the terms in the @pop tab
summary_snps2016 <- summary(snps2016)
#do it by population
summary_Alb <- summary(Albany16)
summary_BF <- summary(BayFarm16)
summary_CB <- summary(CrownBeach16)
summary_Coo <- summary(Coos16)
summary_Hum <- summary(Humboldt16)
summary_Mor <- summary(Morro16)
summary_OS <- summary(OceanShores16)
summary_Par <- summary(Paradise16)
summary_PM <- summary(PtMolate16)
summary_PO <- summary(PtOrient16)
summary_Tom <- summary(Tomales16)

mean(na.omit(summary_Alb$Hobs)) #0.08096898
mean(na.omit(summary_BF$Hobs)) #0.09654527
mean(na.omit(summary_CB$Hobs)) #0.09018269
mean(na.omit(summary_Coo$Hobs)) #0.07572902
mean(na.omit(summary_Hum$Hobs)) #0.07350129
mean(na.omit(summary_Mor$Hobs)) #0.1132373
mean(na.omit(summary_OS$Hobs)) #0.06460843
mean(na.omit(summary_Par$Hobs)) #0.08738069
mean(na.omit(summary_PM$Hobs)) #0.09171222
mean(na.omit(summary_PO$Hobs)) #0.08780595
mean(na.omit(summary_Tom$Hobs)) #0.08274748

#provide summary table in supplement of missing data and total alleles? -- should compare to original # of individuals and total loci in dataset
#Albany: 4 ind, 56383 alleles, 58.12% missing
#Bay Farm: 19 ind, 85470 alleles, 53.10% missing
#Crown Beach: 9 ind, 65479 alleles, 56.44% missing
#Coos: 12 ind, 60532 alleles, 57.59% missing
#Humboldt: 4 ind, 50990 alleles, 59.72% missing
#Morro: 7 ind, 64816 alleles, 57.28% missing
#Ocean Shores: 3 ind, 43744 alleles, 62.30% missing
#Paradise: 7 ind, 64867 alleles, 56.23% missing
#Pt Molate: 14 ind, 71541 alleles, 57.00% missing
#Pt Orient: 6 ind, 60255 alleles, 56.56% missing
#Tomales: 3 ind, 45742 alleles, 61.34% missing

plot(summary_snps2016$Hexp, summary_snps2016$Hobs, pch=20, cex=3, xlim=c(0, 1), ylim=c(0, 1), abline(0,1,lty=2))

#this is just for the whole dataset
#do not use - language defunct most likely 3.29.20
snps2016_pval <- hw.test(snps2016, res="matrix") #still gives warnings even though I changed ploidy - possibly did this wrong?
snps2016_sigPval <- which(snps2016_pval<0.0001,TRUE)
snps2016_pvalAll <- hw.test(snps2016, res="full") #this must take forever to do...
snps2016_completeTest <- mapply(function(i,j) snps2016_pvalAll[[i]][[j]], snps2016_sigPval[,2], snps2016_sigPval[,1], SIMPLIFY=FALSE)

#by population
#Monte Carlo test is not working because it doesn't think they are all diploid - which is not true!!
#B is the number of replicates for the Monte Carlo procedure
#3.29.20 export does not have reliable Monte Carlo results (Pr_exact)
het_Alb_pval <- hw.test(Albany16, B=1000) 
het_Alb_sigPval <- het_Alb_pval[which(het_Alb_pval[,3]<0.05,TRUE),]
write.table(het_Alb_sigPval, "~/Documents/R_Data/ddRADSeq/obsHetSig_Albany_6.15.20.txt", row.names = T, quote = F, sep='\t')
het_BF_pval <- hw.test(BayFarm16, B=1000)
het_BF_sigPval <- het_BF_pval[which(het_BF_pval[,3]<0.05,TRUE),]
write.table(het_BF_sigPval, "~/Documents/R_Data/ddRADSeq/obsHetSig_BayFarm_6.15.20.txt", row.names = T, quote = F, sep='\t')
het_CB_pval <- hw.test(CrownBeach16, B=1000)
het_CB_sigPval <- het_CB_pval[which(het_CB_pval[,3]<0.05,TRUE),]
write.table(het_CB_sigPval, "~/Documents/R_Data/ddRADSeq/obsHetSig_CrownBeach_6.15.20.txt", row.names = T, quote = F, sep='\t')
het_Coo_pval <- hw.test(Coos16, B=1000)
het_Coo_sigPval <- het_Coo_pval[which(het_Coo_pval[,3]<0.05,TRUE),]
write.table(het_Coo_sigPval, "~/Documents/R_Data/ddRADSeq/obsHetSig_Coos_6.15.20.txt", row.names = T, quote = F, sep='\t')
het_Hum_pval <- hw.test(Humboldt16, B=1000)
het_Hum_sigPval <- het_Hum_pval[which(het_Hum_pval[,3]<0.05,TRUE),]
write.table(het_Hum_sigPval, "~/Documents/R_Data/ddRADSeq/obsHetSig_Humboldt_6.15.20.txt", row.names = T, quote = F, sep='\t')
het_Mor_pval <- hw.test(Morro16, B=1000)
het_Mor_sigPval <- het_Mor_pval[which(het_Mor_pval[,3]<0.05,TRUE),]
write.table(het_Mor_sigPval, "~/Documents/R_Data/ddRADSeq/obsHetSig_Morro_6.15.20.txt", row.names = T, quote = F, sep='\t')
het_OS_pval <- hw.test(OceanShores16, B=1000)
het_OS_sigPval <- het_OS_pval[which(het_OS_pval[,3]<0.05,TRUE),]
write.table(het_OS_sigPval, "~/Documents/R_Data/ddRADSeq/obsHetSig_OceanShores_6.15.20.txt", row.names = T, quote = F, sep='\t')
het_Par_pval <- hw.test(Paradise16, B=1000)
het_Par_sigPval <- het_Par_pval[which(het_Par_pval[,3]<0.05,TRUE),]
write.table(het_Par_sigPval, "~/Documents/R_Data/ddRADSeq/obsHetSig_Paradise_6.15.20.txt", row.names = T, quote = F, sep='\t')
het_PM_pval <- hw.test(PtMolate16, B=1000)
het_PM_sigPval <- het_PM_pval[which(het_PM_pval[,3]<0.05,TRUE),]
write.table(het_PM_sigPval, "~/Documents/R_Data/ddRADSeq/obsHetSig_PtMolate_6.15.20.txt", row.names = T, quote = F, sep='\t')
het_PO_pval <- hw.test(PtOrient16, B=1000)
het_PO_sigPval <- het_PO_pval[which(het_PO_pval[,3]<0.05,TRUE),]
write.table(het_PO_sigPval, "~/Documents/R_Data/ddRADSeq/obsHetSig_PtOrient_6.15.20.txt", row.names = T, quote = F, sep='\t')
het_Tom_pval <- hw.test(Tomales16, res="matrix")
het_Tom_sigPval <- het_Tom_pval[which(het_Tom_pval[,3]<0.05,TRUE),]
write.table(het_Tom_sigPval, "~/Documents/R_Data/ddRADSeq/obsHetSig_Tomales_6.15.20.txt", row.names = T, quote = F, sep='\t')

####trees####
#Provesti's distance is specifically for snps
#use the poppr package
library("poppr")
library("ape") # To visualize the tree using the "nj" function
library("magrittr")

#just 2016 dataset
set.seed(999)
snps2016 %>%
  genind2genpop(pop = ~Population) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = provesti.dist)

provDist16 <- provesti.dist(snps2016)
tree16_new <- provDist16 %>%
  bionj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade

tree16 <- provDist16 %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade

#rename tip labels
tree16$tip.label <- as.character(match2016_snpsPOP$ID)
tree16_new$tip.label <- as.character(match2016_snpsPOP$ID)

#how to order tip color based on the order shown in the tree??
#use changeRR21 as the color gradient (red to yellow)
colors_changeRR21

plot(tree16, tip.color = colors_changeRR21)


plot(tree16_new)
