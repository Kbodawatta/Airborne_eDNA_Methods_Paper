###This r script is used to analyse the data from BirT primers for testing effect of sample storage method on accumulation of vertebrate DNA in the study "Advanced airborne eDNA sampling allows robust spatiotemporal characterisation of vertebrate communities"
## Experiment 3_(Tofte Skov)

#Load the packages
library(vegan)
library(ggfortify)
library(devtools)
library(pairwiseAdonis)
library(viridis)
library(ggplot2)
library(ape)
library(dplyr)
library(hagis)
library(ggplot2)
library(lmerTest)
library(emmeans) 

###Analysis of species richness. We are using all samples including samples without detections.
Data <- read.csv2("BirT_Airflow_All.csv", row.names = 1, check.names = FALSE)

#Testing the effect of filter types on species richness using linear mixed models
model1 <- lmer(Richness ~  SampleType +(1|SampleSite) + (1|Date), Data)

summary(model1)
anova(model1)

#Investigating the pairwise differences based on linear mixed models
Pair <- pairs(emmeans(model1,  ~SampleType))
summary(Pair)

#Plotting the richness data with box plots
Data$SampleType <- factor(Data$SampleType, levels = c("12v", "12v_5", "5v", "5v_12", "12vp", "5vp"))

ggplot(Data, aes(x=SampleType, y=Richness)) + geom_boxplot(aes(fill = SampleType), outlier.shape = NA)+  geom_jitter(aes(shape = Date), col = "#b5b1b2", size = 3, width = 0.15) + theme_classic() + scale_fill_manual(values=c("#373252" ,"#cc1d1d","#7a0a0f", "#c70e8c", "#0a0001", "#666262"))


###Species accumulation curve (https://rdrr.io/rforge/vegan/man/specaccum.html). These were done seperately for bestperforming sampling setup and the rest of sampling setups.
Tab5 <- read.csv2("BirT_Others.csv", row.names = 1, check.names = FALSE)
Tab6 <- Tab5[5:40]
sp1 <- specaccum(Tab6)
sp2 <- specaccum(Tab6, "random")
sp2
summary(sp2)
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlim=c(0, 85), ylim=c(0, 42))
boxplot(sp2, col="yellow", add=TRUE, pch="+")

##only for best sampler
Tab3 <- read.csv2("BirT_12v.csv", row.names = 1, check.names = FALSE)
Tab4 <- Tab3[5:39]
sp1 <- specaccum(Tab4)
sp2 <- specaccum(Tab4, "random")
sp2
summary(sp2)
plot(sp1, ci.type="poly", col="red", lwd=2, ci.lty=0, ci.col="lightgray", xlim=c(0, 85), ylim=c(0, 42))
boxplot(sp2, col="lightpink", add=TRUE, pch="+")

### Community level analyses. Here we are only using samples with taxa detections.
#Bring in the datasheet
Tab <- read.csv2("BirT_AirFlow_OTU.csv", row.names = 1, check.names = FALSE)

Tab2 <- Tab[5:45]
Meta <- Tab[1:4]

# making the distance matrix with Raup Crick dissimilarity
dis <- raupcrick(Tab2, null = "r1", nsimul = 999)

#Conducting the PCOA analysis
princoor <- pcoa(dis)

#getting the percentage of variation explained by 1st and 2nd PCOA axis
Axis1.percent <-princoor$values$Relative_eig[[1]] * 100
Axis2.percent <-princoor$values$Relative_eig[[2]] * 100

Axis1.percent
Axis2.percent
princoor.pathotype.data <-data.frame(Sample = as.integer(rownames(princoor$vectors)),X = princoor$vectors[, 1], Y = princoor$vectors[, 2])

#download the coordinates and then add coordinates to the BirT_Airflow_Stat document
write.csv2(princoor.pathotype.data, file = "BirT_Cord_Raup.csv")

#Conduct Permutational multivariate analyses tests (PERMANOVAS) using the adonis2 function. We use the same distance matrix that we used for PCOA.
adonis2(formula = dis~ SampleType + SampleSite, data = Meta, by = "margin", permutations = 10000, strata = Meta$Date)

#If PERMANOVA significant we can use the pairwise.adonis function to test pairwise differences
pairwise.adonis(dis, Meta$SampleType, p.adjust.m = "bonferroni", reduce = NULL,
                perm = 10000)

#Now bring in the stat document with the PCOA coordinates
Data <- read.csv2("BirT_Airflow_Stat.csv", row.names = 1, check.names = FALSE)
Data$SampleType <- factor(Data$SampleType, levels = c("12v", "12v_5", "5v", "5v_12", "12vP"))

#Plot the PCOA plot
ggplot(data = Data, aes(x = Cor_X, y = Cor_Y)) + geom_point(size=7, aes(colour = SampleType, shape = Date))+ scale_color_manual(values=c("#240da3" ,"#9382ed","#7a0a0f", "#e34f56", "#0a0001")) + theme_classic() + scale_shape_manual(values=c(16, 17, 15)) + geom_line(aes(group=SampleSite))
ggplot(data = Data, aes(x = Cor_X, y = Cor_Y)) + geom_point(size=7, aes(colour = SampleType, shape = Date))+ scale_color_manual(values=c("#240da3" ,"#9382ed","#7a0a0f", "#e34f56", "#0a0001")) + theme_classic() + scale_shape_manual(values=c(16, 17, 15)) + stat_ellipse(aes(group = SampleType, col = SampleType), linetype = 2, type = "t", level = 0.80)

