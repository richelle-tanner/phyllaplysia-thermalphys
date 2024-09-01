#R Tanner 12/6/16
#edit for colors and shapes 12.19.20

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

cb_color <- c("#CC79A7", "#D55E00", "#D55E00", "#E69F00", "#E69F00", "#0072B2", "#CC79A7", "#56B4E9", "#009E73", "#0072B2", "#009E73") #alphabetical order

#import data
ddr <- "~/Documents/R_Data/Latitudinal Study/"
envTempDataset <- read.csv(paste(ddr, "LATSF temp vs ARR.csv", sep=""))
CTmaxRRdata <- read.delim("~/Documents/R_Data/ddRADSeq/phenotype_summary_table_6.18.20.txt")

envTempDataset <- envTempDataset[order(envTempDataset$location),]

library(ggplot2)

#plot Avg Temp 2Mcollection vs envARR
lm_eqn_ARR <- function(envTempDataset){
  m <- lm(totalARR ~ Avg.Temp.2Mcollection, envTempDataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
ggplot(data=envTempDataset, aes(x=Avg.Temp.2Mcollection, y=totalARR)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Acclimation Response Ratio") +
  xlab("Average Habitat Temperature (ºC)") +
  #geom_text(x = 15, y = 0.2, label = lm_eqn_ARR(envTempDataset), parse = TRUE) +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not significant
lm_envARR <- lm(totalARR ~ Avg.Temp.2Mcollection, envTempDataset)
summary(lm_envARR)

#plot Avg Temp 2Mcollection vs envCTmax
#get lm equation and r2
lm_eqn <- function(envTempDataset){
  m <- lm(relevant_envCTmax ~ Avg.Temp.2Mcollection, envTempDataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

#split for use of two different shapes in graph
envTempDataset_LAT <- envTempDataset[c(3,5,6,7,11),]
envTempDataset_SF <- envTempDataset[-c(3,5,6,7,11),]

ggplot(data=envTempDataset_LAT, aes(x=Avg.Temp.2Mcollection, y=relevant_envCTmax)) +
  geom_point(size=10, shape=24, aes(fill=location)) +
  geom_point(data=envTempDataset_SF, size=10, shape=23, aes(fill=location)) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Critical Thermal Maximum (ºC)") +
  xlab("Average Habitat Temperature (ºC)") +
  geom_text(x = 14, y = 32.5, label = lm_eqn(envTempDataset), parse = TRUE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_fill_manual(values=cb_color)
#ANOVA - significant*****p<0.05
lm_envCTmax <- lm(relevant_envCTmax ~ Avg.Temp.2Mcollection, envTempDataset)
summary(lm_envCTmax)

#plot CC.Projection  vs envARR
ggplot(data=envTempDataset, aes(x=CC.Projection, y=relevant_envARR)) +
  geom_point(size=5) +
  geom_smooth(method='lm')
#ANOVA - not significant
lm_CCenvARR <- lm(relevant_envARR ~ CC.Projection, envTempDataset)
summary(lm_CCenvARR)

#plot CC.Projection vs. envCTmax  
ggplot(data=envTempDataset, aes(x=CC.Projection, y=relevant_envCTmax)) +
  geom_point(size=5) +
  geom_smooth(method='lm')
#ANOVA - significant******p<0.05
lm_CCenvCTmax <- lm(relevant_envCTmax ~ CC.Projection, envTempDataset)
summary(lm_CCenvCTmax)

#plot Daily Variation vs envCTmax
lm_eqn_CTmax <- function(envTempDataset){
  m <- lm(X13CTmax ~ Daily.Variation, envTempDataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
ggplot(data=envTempDataset_LAT, aes(x=Daily.Variation, y=X13CTmax)) +
  geom_point(size=10, shape=24, aes(fill=location)) +
  geom_point(data=envTempDataset_SF, size=10, shape=23, aes(fill=location)) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Critical Thermal Maximum (ºC)") +
  xlab("Average Daily Variation in Temperature (ºC)") +
  geom_text(x = 6, y = 31, label = lm_eqn_CTmax(envTempDataset), parse = TRUE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_fill_manual(values=cb_color)
  
#ANOVA -  significant
lm_DVenvCTmax <- lm(X13CTmax ~ Daily.Variation, envTempDataset)
summary(lm_DVenvCTmax)

#plot Daily Variation vs ARR
lm_eqn_ARR <- function(envTempDataset){
  m <- lm(totalARR ~ Daily.Variation, envTempDataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
ggplot(data=envTempDataset, aes(x=Daily.Variation, y=totalARR)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Acclimation Response Ratio") +
  xlab("Average Daily Variation in Temperature (ºC)") +
  #geom_text(x = 4, y = 0.8, label = lm_eqn_ARR(envTempDataset), parse = TRUE) +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA
lm_DVenvARR <- lm(totalARR ~ Daily.Variation, envTempDataset)
summary(lm_DVenvARR)



#ARR vs latitude
lm_eqn_ARR <- function(envTempDataset){
  m <- lm(totalARR ~ Latitude, envTempDataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
ggplot(data=envTempDataset, aes(x=Latitude, y=totalARR)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Acclimation Response Ratio") +
  xlab("Latitude") +
  geom_text(x = 42, y = 0.6, label = lm_eqn_ARR(envTempDataset), parse = TRUE) +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))

ggplot(data=envTempDataset, aes(x=Latitude, y=X13CTmax)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Critical Thermal Maximum") +
  xlab("Latitude") +
  #geom_text(x = 15, y = 0.2, label = lm_eqn_ARR(envTempDataset), parse = TRUE) +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))


#temp13 and temp13rr come from CTmax vs resp rate.R
#change13 <- merge(envTempDataset, temp13, by = "location")
#RR13 <- merge(envTempDataset, temp13rr, by = "location")

#now use the new one CTmaxRRdata (imported at top)
colnames(CTmaxRRdata)[7] <- "location"
RR <- merge(envTempDataset, CTmaxRRdata, by = "location")

ggplot(data=RR, aes(x=Avg.Temp.2Mcollection, y=preHS_RR_13avg)) +
  geom_point(size=15, aes(color=location)) +
  #geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Average Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_respRateAvgTemp13 <- lm(preHS_RR_13avg ~ Avg.Temp.2Mcollection, RR)
summary(lm_respRateAvgTemp13)

ggplot(data=RR, aes(x=Avg.Temp.2Mcollection, y=preHS_RR_17avg)) +
  geom_point(size=5) +
  #geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Average Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_respRateAvgTemp17 <- lm(preHS_RR_17avg ~ Avg.Temp.2Mcollection, RR)
summary(lm_respRateAvgTemp17)

ggplot(data=RR, aes(x=Avg.Temp.2Mcollection, y=preHS_RR_21avg)) +
  geom_point(size=5) +
  #geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Average Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_respRateAvgTemp21 <- lm(preHS_RR_21avg ~ Avg.Temp.2Mcollection, RR)
summary(lm_respRateAvgTemp21)



ggplot(data=change13, aes(x=Avg.Temp.2Mcollection, y=change)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Change in Respiration Rate") +
  xlab("Average Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_changeAvgTemp <- lm(change ~ Avg.Temp.2Mcollection, change13)
summary(lm_changeAvgTemp)



ggplot(data=change13, aes(x=Daily.Variation, y=change)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Change in Respiration Rate") +
  xlab("Average Daily Variation in Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_changeDV <- lm(change ~ Daily.Variation, change13)
summary(lm_changeDV)

ggplot(data=RR13, aes(x=Daily.Variation, y=respRate)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Average Daily Variation in Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_respRateDV <- lm(respRate ~ Daily.Variation, RR13)
summary(lm_respRateDV)

ggplot(data=change13, aes(x=Latitude, y=change)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Change in Respiration Rate") +
  xlab("Latitude") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_changeLat <- lm(change ~ Latitude, change13)
summary(lm_changeLat)

ggplot(data=RR13, aes(x=Latitude, y=respRate)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Latitude") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_respRateLat <- lm(respRate ~ Latitude, RR13)
summary(lm_respRateLat)

#what about 21ºC respiration rates?
#still uses CTmax vs resp rates.R
change21 <- merge(envTempDataset, temp21, by = "location")
RR21 <- merge(envTempDataset, temp21rr, by = "location")

ggplot(data=change21, aes(x=Avg.Temp.2Mcollection, y=change)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Change in Respiration Rate") +
  xlab("Average Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - sig
lm_changeAvgTemp <- lm(change ~ Avg.Temp.2Mcollection, change21)
summary(lm_changeAvgTemp)

ggplot(data=RR21, aes(x=Avg.Temp.2Mcollection, y=respRate)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Average Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_respRateAvgTemp <- lm(respRate ~ Avg.Temp.2Mcollection, RR21)
summary(lm_respRateAvgTemp)

ggplot(data=change21, aes(x=Daily.Variation, y=change)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Change in Respiration Rate") +
  xlab("Average Daily Variation in Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_changeDV <- lm(change ~ Daily.Variation, change21)
summary(lm_changeDV)

ggplot(data=RR21, aes(x=Daily.Variation, y=respRate)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Average Daily Variation in Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_respRateDV <- lm(respRate ~ Daily.Variation, RR21)
summary(lm_respRateDV)

ggplot(data=change21, aes(x=Latitude, y=change)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Change in Respiration Rate") +
  xlab("Latitude") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_changeLat <- lm(change ~ Latitude, change21)
summary(lm_changeLat)

ggplot(data=RR21, aes(x=Latitude, y=respRate)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Latitude") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_respRateLat <- lm(respRate ~ Latitude, RR21)
summary(lm_respRateLat)

#what about 17ºC respiration rates?
#still uses CTmax vs resp rates.R
change17 <- merge(envTempDataset, temp17, by = "location")
RR17 <- merge(envTempDataset, temp17rr, by = "location")

ggplot(data=change17, aes(x=Avg.Temp.2Mcollection, y=change)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Change in Respiration Rate") +
  xlab("Average Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - sig
lm_changeAvgTemp <- lm(change ~ Avg.Temp.2Mcollection, change17)
summary(lm_changeAvgTemp)

ggplot(data=RR17, aes(x=Avg.Temp.2Mcollection, y=respRate)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Average Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_respRateAvgTemp <- lm(respRate ~ Avg.Temp.2Mcollection, RR17)
summary(lm_respRateAvgTemp)

ggplot(data=change17, aes(x=Daily.Variation, y=change)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Change in Respiration Rate") +
  xlab("Average Daily Variation in Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_changeDV <- lm(change ~ Daily.Variation, change17)
summary(lm_changeDV)

ggplot(data=RR17, aes(x=Daily.Variation, y=respRate)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Average Daily Variation in Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_respRateDV <- lm(respRate ~ Daily.Variation, RR17)
summary(lm_respRateDV)

ggplot(data=change17, aes(x=Latitude, y=change)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Change in Respiration Rate") +
  xlab("Latitude") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_changeLat <- lm(change ~ Latitude, change17)
summary(lm_changeLat)

ggplot(data=RR17, aes(x=Latitude, y=respRate)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Latitude") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_respRateLat <- lm(respRate ~ Latitude, RR17)
summary(lm_respRateLat)




#absolute maximum temperature 
ggplot(data=envTempDataset, aes(x=abs_maxTemp, y=totalARR)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Acclimation Response Ratio") +
  xlab("Maximum Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_maxARR <- lm(totalARR ~ abs_maxTemp, envTempDataset)
summary(lm_maxARR)

ggplot(data=envTempDataset, aes(x=abs_maxTemp, y=X13CTmax)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Critical Thermal Maximum") +
  xlab("Maximum Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_maxCTmax13 <- lm(X13CTmax ~ abs_maxTemp, envTempDataset)
summary(lm_maxCTmax13)

ggplot(data=change13, aes(x=abs_maxTemp, y=change)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Change in Respiration Rate") +
  xlab("Maximum Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_maxChange <- lm(change ~ abs_maxTemp, change13)
summary(lm_maxChange)

ggplot(data=change21, aes(x=abs_maxTemp, y=change)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Change in Respiration Rate") +
  xlab("Maximum Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_maxChange <- lm(change ~ abs_maxTemp, change21)
summary(lm_maxChange)

ggplot(data=RR13, aes(x=abs_maxTemp, y=respRate)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Maximum Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_maxRespRate <- lm(respRate ~ abs_maxTemp, RR13)
summary(lm_maxRespRate)

ggplot(data=RR21, aes(x=abs_maxTemp, y=respRate)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Resting Respiration Rate") +
  xlab("Maximum Habitat Temperature (ºC)") +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - not sig
lm_maxRespRate <- lm(respRate ~ abs_maxTemp, RR21)
summary(lm_maxRespRate)










#plot Avg Temp 2Mcollection vs envCTmax
#get lm equation and r2
lm_eqn <- function(envTempDataset){
  m <- lm(X13CTmax ~ Avg.Temp.2Mcollection, envTempDataset);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
ggplot(data=envTempDataset, aes(x=Avg.Temp.2Mcollection, y=X13CTmax)) +
  geom_point(size=5) +
  geom_smooth(method='lm', se=FALSE) +
  theme_bw() +
  ylab("Critical Thermal Maximum (ºC)") +
  xlab("Average Habitat Temperature (ºC)") +
  geom_text(x = 14, y = 30, label = lm_eqn(envTempDataset), parse = TRUE) +
  theme(axis.text.x = element_text(size=15, face="bold"),
        axis.text.y = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))
#ANOVA - significant*****p<0.05
lm_envCTmax <- lm(X13CTmax ~ Avg.Temp.2Mcollection, envTempDataset)
summary(lm_envCTmax)


####edit 6.15.20 for new figure####
#add in rows that just have hab temp if they don't have change in RR at 21
add_rows <- envTempDataset[c(7,10,11),]
add_rows$temp <- 21
add_rows$N <- 0
add_rows$change <- 0
add_rows$sd <- 0
add_rows$se <- 0
add_rows$ci <- 0
change21 <- rbind(change21, add_rows)
write.table(change21, "~/Documents/R_Data/ddRADSeq/habTemp_changeRR21_6.15.20.txt", quote=F, sep="\t", row.names=F)


#plot the points so we can put them on the map in illustrator
ggplot(change21, aes(x=location, y=Latitude)) +
  geom_point(aes(size=Avg.Temp.2Mcollection, color=change)) +
  scale_color_gradient(low="blue", high="red")
