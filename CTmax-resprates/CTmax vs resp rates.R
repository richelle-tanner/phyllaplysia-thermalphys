#plot CTmax and change in resp rate for each individual against each other, color by location
#major edit June 2020 for new figures and tables - use this script for everything phenotype-related (CTmax, resp rates, ARR)

#color codes for each location:
#Color Blind Friendly Palette: cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#Albany: #CC79A7
#Bay Farm: #D55E00
#Coos Bay: #D55E00
#Crown Beach: #E69F00
#Humboldt Bay: #E69F00
#Morro Bay: #0072B2
#Ocean Shores: #CC79A7
#Paradise: #56B4E9
#Point Molate: #009E73
#Point Orient: #0072B2
#Tomales Bay: #009E73

#packages
library(ggplot2)

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#import data
ddr <- "~/Documents/R_Data/Latitudinal Study/"
respRates <- read.csv(paste(ddr, "respiration rates.csv", sep=""))
CTmax <- read.csv(paste(ddr, "CTmax_allPTay.csv", sep=""))

#reformat columns
CTmax$vial <- sub("^[^_]*", "", CTmax$Full_ID)
CTmax$vial <- gsub("_", '', CTmax$vial)
CTmax$group <- paste(CTmax$Group, CTmax$Acclim_temp, "_", CTmax$vial, sep="")
#just a dataframe with CTmax and group EDIT 2.23.20 added ARR and weight so we can have a comprehensive data frame for export to adegenet
CTmax1 <- data.frame(treatment=CTmax$group, location=CTmax$Location, temp=CTmax$Acclim_temp, 
                     CTmax=CTmax$Ctmax_degC, ARR=CTmax$ARR_overall, weight_mg=CTmax$Weight_mg)

#subtract noHS from HS
#make new individual ID column from vial and treatment
respRates$ID <- paste(respRates$group, "_", respRates$vial, sep="")
#subset HS and noHS dataframes
HS <- subset(respRates, HS == "HS")
noHS <- subset(respRates, HS == "noHS")
changeRR <- data.frame(treatment=HS$ID, group=HS$group, vial=HS$vial, respRate_noHS = noHS$mgO2hrg, 
                       respRate_HS = HS$mgO2hrg, change = HS$mgO2hrg - noHS$mgO2hrg)
combinedData <- merge(changeRR, CTmax1, by = "treatment")

#export the combined data frame to use with genetic info EDIT 2.23.20
write.table(combinedData, "~/Documents/R_Data/ddRADSeq/combined_respRateCTmax.txt", quote = F, sep='\t', row.names = F)

#summarize respiration rates
summaryChange <- summarySE(changeRR, measurevar = "change", groupvars = "group")
summaryChangeLocation <- summarySE(combinedData, measurevar = "change", groupvars = "location")

# #summarize data
summaryChange <- summarySE(combinedData, measurevar = "change", groupvars = c("temp", "location"))
summaryRR <- summarySE(combinedData, measurevar = "respRate", groupvars = c("temp", "location"))

temp13 <- subset(summaryChange, temp == 13)
temp13rr <- subset(summaryRR, temp == 13)
temp17 <- subset(summaryChange, temp == 17)
temp17rr <- subset(summaryRR, temp == 17)
temp21 <- subset(summaryChange, temp == 21)
temp21rr <- subset(summaryRR, temp == 21)

CTmaxVrr <- aov(CTmax ~ respRate_noHS + location, data=combinedData)
summary(CTmaxVrr)

CTmaxVrrHS <- aov(CTmax ~ respRate_HS + location, data=combinedData)
summary(CTmaxVrrHS)

#Albany: #CC79A7
#Bay Farm: #D55E00
#Coos Bay: #D55E00
#Crown Beach: #E69F00
#Humboldt Bay: #E69F00
#Morro Bay: #0072B2
#Ocean Shores: #CC79A7
#Paradise: #56B4E9
#Point Molate: #009E73
#Point Orient: #0072B2
#Tomales Bay: #009E73

ggplot(combinedData, aes(x=respRate_HS, y=CTmax, group=location, fill=location, shape=location)) +
  geom_point(aes(color=location, size=5)) +
  ylab("Critical thermal limits (°C)") +
  xlab("Respiration Rate After HS (mg O2/hr/g)") +
  theme_bw() +
  scale_color_manual(values=c("#CC79A7", "#D55E00", "#D55E00", "#E69F00", "#E69F00", "#0072B2",
                              "#CC79A7", "#56B4E9", "#009E73", "#0072B2", "#009E73")) +
  scale_fill_manual(values=c("#CC79A7", "#D55E00", "#D55E00", "#E69F00", "#E69F00", "#0072B2",
                              "#CC79A7", "#56B4E9", "#009E73", "#0072B2", "#009E73")) +
  scale_shape_manual(values=c(23,23,17,23,17,17,17,23,23,23,17))

ggplot(combinedData, aes(x=respRate_noHS, y=CTmax, group=location, fill=location, shape=location)) +
  geom_point(aes(color=location, size=5)) +
  ylab("Critical thermal limits (°C)") +
  xlab("Respiration Rate Before HS (mg O2/hr/g)") +
  theme_bw() +
  scale_color_manual(values=c("#CC79A7", "#D55E00", "#D55E00", "#E69F00", "#E69F00", "#0072B2",
                              "#CC79A7", "#56B4E9", "#009E73", "#0072B2", "#009E73")) +
  scale_fill_manual(values=c("#CC79A7", "#D55E00", "#D55E00", "#E69F00", "#E69F00", "#0072B2",
                             "#CC79A7", "#56B4E9", "#009E73", "#0072B2", "#009E73")) +
  scale_shape_manual(values=c(23,23,17,23,17,17,17,23,23,23,17))

ggplot(combinedData, aes(x=change, y=CTmax, group=location, fill=location, shape=location)) +
  geom_point(aes(color=location, size=5)) +
  ylab("Critical thermal limits (°C)") +
  xlab("Change in Respiration Rate Pre/Post HS (mg O2/hr/g)") +
  theme_bw() +
  scale_color_manual(values=c("#CC79A7", "#D55E00", "#D55E00", "#E69F00", "#E69F00", "#0072B2",
                              "#CC79A7", "#56B4E9", "#009E73", "#0072B2", "#009E73")) +
  scale_fill_manual(values=c("#CC79A7", "#D55E00", "#D55E00", "#E69F00", "#E69F00", "#0072B2",
                             "#CC79A7", "#56B4E9", "#009E73", "#0072B2", "#009E73")) +
  scale_shape_manual(values=c(23,23,17,23,17,17,17,23,23,23,17))

####new visualizations June 2020####

#recalculate ARR 
#use max-min idea (from J 6.15.20)
combinedData$ARR <- NULL #get rid of old calculation
CTmax1$ARR <- NULL
#collect summarized data
CTmax_summary <- matrix(nrow = 11, ncol = 6)
ARR <- vector()
#first summarize by location
pop_names <- as.character(unique(combinedData$location))
for (i in 1:length(pop_names)) {
  #then summarize by acclimation
  temp_pop <- subset(CTmax1, location == pop_names[i])
  treat <- as.character(unique(CTmax1$temp))
  #then take average and SD of max and min
  #then subtract and populate SD error
  if (length(treat) == "3") {
    temp_pop13 <- subset(temp_pop, temp == 13)
    temp_pop17 <- subset(temp_pop, temp == 17)
    temp_pop21 <- subset(temp_pop, temp == 21)
    avg13 <- mean(temp_pop13$CTmax)
    sd13 <- sd(temp_pop13$CTmax)
    avg17 <- mean(temp_pop17$CTmax)
    sd17 <- sd(temp_pop17$CTmax)
    avg21 <- mean(temp_pop21$CTmax)
    sd21 <- sd(temp_pop21$CTmax)
    CTmax_summary[i,] <- c(avg13, sd13, avg17, sd17, avg21, sd21)
    ARR[i] <- (avg21 - avg13)/8
  }
  if (length(treat) == "2") {
    temp_pop13 <- subset(temp_pop, temp == 13)
    temp_pop17 <- subset(temp_pop, temp == 17)
    avg13 <- mean(temp_pop13$CTmax)
    sd13 <- sd(temp_pop13$CTmax)
    avg17 <- mean(temp_pop17$CTmax)
    sd17 <- sd(temp_pop17$CTmax)
    CTmax_summary[i,] <- c(avg13, sd13, avg17, sd17, "NA", "NA")
    ARR[i] <- (avg17 - avg13)/4
  }
  if (length(treat) == "1") {
    temp_pop13 <- subset(temp_pop, temp == 13)
    avg13 <- mean(temp_pop13$CTmax)
    sd13 <- sd(temp_pop13$CTmax)
    CTmax_summary[i,] <- c(avg13, sd13, "NA", "NA", "NA", "NA")
    ARR[i] <- "NA"
  }
}

CTmax_summary <- data.frame(CTmax_summary)
colnames(CTmax_summary) <- c("CTmax_13avg", "CTmax_13sd", "CTmax_17avg", "CTmax_17sd", "CTmax_21avg", "CTmax_21sd")
CTmax_summary$population <- pop_names

changeRR$group <- as.character(changeRR$group)
changeRR$location <- substr(changeRR$group, 0, nchar(changeRR$group)-2)
changeRR$temp <- as.numeric(substr(changeRR$group, nchar(changeRR$group)-1, nchar(changeRR$group)))
RR_summary <- matrix(nrow = 11, ncol = 12)
RR_plasticityPre <- vector()
RR_plasticityPost <- vector()
pop_namesRR <- as.character(unique(changeRR$location))
#first summarize by location
for (i in 1:length(pop_namesRR)) {
  #then summarize by acclimation
  temp_pop <- subset(changeRR, location == pop_namesRR[i])
  treat <- as.character(unique(temp_pop$temp))
  #then take average and SD of max and min
  #then subtract and populate SD error
  if (length(treat) == "3") {
    temp_pop13 <- subset(temp_pop, temp == 13)
    temp_pop17 <- subset(temp_pop, temp == 17)
    temp_pop21 <- subset(temp_pop, temp == 21)
    avg13pre <- mean(temp_pop13$respRate_noHS)
    sd13pre <- sd(temp_pop13$respRate_noHS)
    avg17pre <- mean(temp_pop17$respRate_noHS)
    sd17pre <- sd(temp_pop17$respRate_noHS)
    avg21pre <- mean(temp_pop21$respRate_noHS)
    sd21pre <- sd(temp_pop21$respRate_noHS)
    avg13post <- mean(temp_pop13$respRate_HS)
    sd13post <- sd(temp_pop13$respRate_HS)
    avg17post <- mean(temp_pop17$respRate_HS)
    sd17post <- sd(temp_pop17$respRate_HS)
    avg21post <- mean(temp_pop21$respRate_HS)
    sd21post <- sd(temp_pop21$respRate_HS)
    RR_summary[i,] <- c(avg13pre, sd13pre, avg17pre, sd17pre, avg21pre, sd21pre, avg13post, sd13post, avg17post, sd17post, avg21post, sd21post)
    RR_plasticityPre[i] <- (avg21pre - avg13pre)/8
    RR_plasticityPost[i] <- (avg21post - avg13post)/8
  }
  if (length(treat) == "2") {
    temp_pop13 <- subset(temp_pop, temp == 13)
    temp_pop17 <- subset(temp_pop, temp == 17)
    avg13pre <- mean(temp_pop13$respRate_noHS)
    sd13pre <- sd(temp_pop13$respRate_noHS)
    avg17pre <- mean(temp_pop17$respRate_noHS)
    sd17pre <- sd(temp_pop17$respRate_noHS)
    avg13post <- mean(temp_pop13$respRate_HS)
    sd13post <- sd(temp_pop13$respRate_HS)
    avg17post <- mean(temp_pop17$respRate_HS)
    sd17post <- sd(temp_pop17$respRate_HS)
    RR_summary[i,] <- c(avg13pre, sd13pre, avg17pre, sd17pre, "NA", "NA", avg13post, sd13post, avg17post, sd17post, "NA", "NA")
    RR_plasticityPre[i] <- (avg17pre - avg13pre)/4
    RR_plasticityPost[i] <- (avg17post - avg13post)/4
  }
  if (length(treat) == "1") {
    if (treat[1] == 13) {
      temp_pop13 <- subset(temp_pop, temp == 13)
      avg13pre <- mean(temp_pop13$respRate_noHS)
      sd13pre <- sd(temp_pop13$respRate_noHS)
      avg13post <- mean(temp_pop13$respRate_HS)
      sd13post <- sd(temp_pop13$respRate_HS)
      RR_summary[i,] <- c(avg13pre, sd13pre, "NA", "NA", "NA", "NA", avg13post, sd13post, "NA", "NA", "NA", "NA")
      RR_plasticityPre[i] <- "NA"
      RR_plasticityPost[i] <- "NA"
    }
    if (treat[1] == 17) {
      temp_pop17 <- subset(temp_pop, temp == 17)
      avg17pre <- mean(temp_pop17$respRate_noHS)
      sd17pre <- sd(temp_pop17$respRate_noHS)
      avg17post <- mean(temp_pop17$respRate_HS)
      sd17post <- sd(temp_pop17$respRate_HS)
      RR_summary[i,] <- c("NA", "NA", avg17pre, sd17pre, "NA", "NA", "NA", "NA", avg17post, sd17post, "NA", "NA")
      RR_plasticityPre[i] <- "NA"
      RR_plasticityPost[i] <- "NA"
    }
  }
}

RR_summary <- data.frame(RR_summary)
colnames(RR_summary) <- c("preHS_RR_13avg", "preHS_RR_13sd", "preHS_RR_17avg", "preHS_RR_17sd", "preHS_RR_21avg", "preHS_RR_21sd",
                          "postHS_RR_13avg", "postHS_RR_13sd", "postHS_RR_17avg", "postHS_RR_17sd", "postHS_RR_21avg", "postHS_RR_21sd")
RR_summary$population <- pop_namesRR
RR_summary$population[7] <- "PtOrient"

#comprehensive table
#columns: population, acclimation temperature, CTmax, pre-HS resp rate, post-HS resp rate, ARR
RR_summary <- RR_summary[order(RR_summary$population),]
CTmax_summary <- CTmax_summary[order(CTmax_summary$population),]
population <- RR_summary$population
cb_color <- c("#CC79A7", "#D55E00", "#D55E00", "#E69F00", "#E69F00", "#0072B2", "#CC79A7", "#56B4E9", "#009E73", "#0072B2", "#009E73") #alphabetical order
RR_summary$population <- NULL

all_summary <- cbind(CTmax_summary, RR_summary)
write.table(all_summary, "~/Documents/R_Data/ddRADSeq/phenotype_summary_table_6.18.20.txt", quote=F, row.names=F, sep='\t')

#figure with plasticity vs trait
#ARR[1] <- 0.06 had to do this manually because of some weird code issues
#ARR[8] <- 0.802083333 same here
#NONE OF THESE ARE SIGNIFICANT - can't talk about beneficial acclimation hypothesis
#plasticity against 13C acclimation
plot_plasticity13CTmax <- data.frame(CTmax=CTmax_summary$CTmax_13avg, population=CTmax_summary$population, ARR=ARR)
summary(lm(ARR ~ CTmax, plot_plasticity13CTmax))
ggplot(plot_plasticity13CTmax, aes(x=CTmax, y=ARR)) +
  geom_point() +
  theme_bw()
#plasticity against 17C acclimation
plot_plasticity17CTmax <- data.frame(CTmax=CTmax_summary$CTmax_17avg, population=CTmax_summary$population, ARR=ARR)
summary(lm(ARR ~ CTmax, plot_plasticity17CTmax))
ggplot(plot_plasticity17CTmax, aes(x=CTmax, y=ARR)) +
  geom_point() +
  theme_bw()
#plasticity against 21C acclimation
plot_plasticity21CTmax <- data.frame(CTmax=CTmax_summary$CTmax_21avg, population=CTmax_summary$population, ARR=ARR)
summary(lm(ARR ~ CTmax, plot_plasticity21CTmax))
#NEW 12.10.20 separate out coastal and SF bay so points can be different shapes
plot_plasticity21CTmax_LAT <- plot_plasticity21CTmax[c(3,5,6,7,11),]
plot_plasticity21CTmax_SF <- plot_plasticity21CTmax[-c(3,5,6,7,11),]
ggplot(plot_plasticity21CTmax_LAT, aes(x=CTmax, y=ARR)) +
  geom_point(size=8, shape=24, aes(fill=population)) +
  geom_point(data=plot_plasticity21CTmax_SF, size=8, shape=23, aes(fill=population)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_fill_manual(values=cb_color)
# Call:
#   lm(formula = ARR ~ CTmax, data = plot_plasticity21CTmax)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.20960 -0.12609  0.00789  0.07304  0.24244 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.44782    0.96762  -0.463    0.660
# CTmax        0.02139    0.03102   0.689    0.516
# 
# Residual standard error: 0.165 on 6 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.07338,	Adjusted R-squared:  -0.08106 
# F-statistic: 0.4751 on 1 and 6 DF,  p-value: 0.5164

#21C is significant!
#plasticity against 13C acclimation (pre HS)
plot_plasticity13RRpre <- data.frame(RR=as.numeric(as.character(RR_summary$preHS_RR_13avg)), population=population, plasticity=as.numeric(RR_plasticityPre))
summary(lm(plasticity ~ RR, plot_plasticity13RRpre))
ggplot(plot_plasticity13RRpre, aes(x=RR, y=plasticity)) +
  geom_point() +
  theme_bw()
#plasticity against 17C acclimation (pre HS)
plot_plasticity17RRpre <- data.frame(RR=as.numeric(as.character(RR_summary$preHS_RR_17avg)), population=population, plasticity=as.numeric(RR_plasticityPre))
summary(lm(plasticity ~ RR, plot_plasticity17RRpre))
ggplot(plot_plasticity17RRpre, aes(x=RR, y=plasticity)) +
  geom_point() +
  theme_bw()
#plasticity against 21C acclimation (pre HS)
plot_plasticity21RRpre <- data.frame(RR=as.numeric(as.character(RR_summary$preHS_RR_21avg)), population=population, plasticity=as.numeric(RR_plasticityPre))
summary(lm(plasticity ~ RR, plot_plasticity21RRpre))
lm_eqn <- function(data){
  m <- lm(data$plasticity ~ data$RR, data);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

plot_plasticity21RRpre_LAT <- plot_plasticity21RRpre[c(3,5,6,7,11),]
plot_plasticity21RRpre_SF <- plot_plasticity21RRpre[-c(3,5,6,7,11),]

ggplot(plot_plasticity21RRpre_LAT, aes(x=RR, y=plasticity)) +
  geom_point(size=8, shape=24, aes(fill=population)) +
  geom_point(data=plot_plasticity21RRpre_SF, size=8, shape=23, aes(fill=population)) +
  theme_bw() +
  geom_text(x = 2, y = 0.02, label = lm_eqn(plot_plasticity21RRpre), parse = TRUE) +
  geom_smooth(method='lm', se=FALSE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_fill_manual(values=cb_color)
  

# Call:
#   lm(formula = plasticity ~ RR, data = plot_plasticity21RRpre)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.081978 -0.024508  0.007611  0.027385  0.065305 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept) -0.05689    0.04471  -1.272   0.2503  
# RR           0.07385    0.02670   2.765   0.0326 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05097 on 6 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.5604,	Adjusted R-squared:  0.4871 
# F-statistic: 7.648 on 1 and 6 DF,  p-value: 0.03262

#21C is significant but heavily driven by one point
#plasticity against 13C acclimation (post HS)
plot_plasticity13RRpost <- data.frame(RR=as.numeric(as.character(RR_summary$postHS_RR_13avg)), population=population, plasticity=as.numeric(RR_plasticityPost))
summary(lm(plasticity ~ RR, plot_plasticity13RRpost))
ggplot(plot_plasticity13RRpost, aes(x=RR, y=plasticity)) +
  geom_point() +
  theme_bw()
#plasticity against 17C acclimation (post HS)
plot_plasticity17RRpost <- data.frame(RR=as.numeric(as.character(RR_summary$postHS_RR_17avg)), population=population, plasticity=as.numeric(RR_plasticityPost))
summary(lm(plasticity ~ RR, plot_plasticity17RRpost))
ggplot(plot_plasticity17RRpost, aes(x=RR, y=plasticity)) +
  geom_point() +
  theme_bw()
#plasticity against 21C acclimation (post HS)
RR_summary$postHS_RR_21avg <- as.numeric(as.character(RR_summary$postHS_RR_21avg))
plot_plasticity21RRpost <- data.frame(RR=RR_summary$postHS_RR_21avg, population=population, plasticity=as.numeric(RR_plasticityPost))
summary(lm(plasticity ~ RR, plot_plasticity21RRpost))
ggplot(plot_plasticity21RRpost, aes(x=RR, y=plasticity)) +
  geom_point(size=8, shape=23, stroke=2, aes(fill=population)) +
  theme_bw() +
  geom_text(x = 7, y = 0.05, label = lm_eqn(plot_plasticity21RRpost), parse = TRUE) +
  geom_smooth(method='lm', se=FALSE)
# Call:
#   lm(formula = plasticity ~ RR, data = plot_plasticity21RRpost)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.46064  0.04335  0.06178  0.07942  0.10591 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.10730    0.08474  -1.266 0.252341    
# RR           0.11907    0.01477   8.059 0.000195 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2031 on 6 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.9154,	Adjusted R-squared:  0.9013 
# F-statistic: 64.94 on 1 and 6 DF,  p-value: 0.0001954


####plot [raw] CTmax values by population####
library(reshape2)
CTmax_avgs <- CTmax_summary[,c(1,3,5,7)]
CTmax_sds <- CTmax_summary[,c(2,4,6,7)]
plot_CTmax_avgs <- melt(CTmax_avgs)
plot_CTmax_sds <- melt(CTmax_sds)
plot_CTmax_summary <- cbind(plot_CTmax_avgs, plot_CTmax_sds)
plot_CTmax_summary$acclimation <- substr(plot_CTmax_summary$variable, 7, 8)
plot_CTmax_summary <- plot_CTmax_summary[,c(1,3,6,7)]
colnames(plot_CTmax_summary) <- c("population", "mean", "SD", "acclimation")
plot_CTmax_summary$SD[is.na(plot_CTmax_summary$SD)] <- 0

#separate out SF bay and coastal populations
plot_CTmax_sum_SF <- plot_CTmax_summary[c(1,2,4,8,9,10,12,13,15,19,20,21,23,24,26,30,31,32),]
plot_CTmax_sum_LAT <- plot_CTmax_summary[-c(1,2,4,8,9,10,12,13,15,19,20,21,23,24,26,30,31,32),]

library(ggplot2)

cb_color_LAT <- c("#D55E00", "#E69F00", "#0072B2", "#CC79A7", "#009E73")
cb_color_SF <- c("#CC79A7", "#D55E00", "#E69F00", "#56B4E9", "#009E73", "#0072B2")

pd <- position_dodge(width = 0.4)
ggplot(plot_CTmax_sum_SF, aes(x=acclimation, y=mean, group=population)) +
  geom_errorbar(aes(x=acclimation, ymin=mean-SD, ymax=mean+SD), width=1, size=1, position=pd) +
  geom_line(position=pd, size=1) +
  geom_point(size=8, aes(group=population, fill=population), shape=23, stroke=1, position=pd) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_fill_manual(values=cb_color_SF)

pd <- position_dodge(width = 0.4)
ggplot(plot_CTmax_sum_LAT, aes(x=acclimation, y=mean, group=population)) +
  geom_errorbar(aes(x=acclimation, ymin=mean-SD, ymax=mean+SD), width=1, size=1, position=pd) +
  geom_line(position=pd, size=1) +
  geom_point(size=8, aes(group=population, fill=population), shape=24, stroke=1, position=pd) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_fill_manual(values=cb_color_LAT)
