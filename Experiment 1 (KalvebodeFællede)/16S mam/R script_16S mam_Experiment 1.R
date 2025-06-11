###This r script is used to analyse the data from 16S mam primers for testing effect of filters on accumulation of vertebrate DNA in the study "Advanced airborne eDNA sampling allows robust spatiotemporal characterisation of vertebrate communities"
## Experiment 1_(KalvebodeFællede)

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
Rich <- read.csv2("16S_Stat_all_sam.csv", row.names = 1, check.names = FALSE)

#Testing the effect of filter types on species richness using linear mixed models
model1 <- lmer(Richness ~  Filter_Type + (1|Date) , Rich)
anova(model1)

#Investigating the pairwise differences based on linear mixed models
Pair <- pairs(emmeans(model1,  ~Filter_Type))
summary(Pair)

#Plotting the richness data with box plots
Rich$Filter_Type <- factor(Rich$Filter_Type, levels = c("F8", "F7_FB", "F7_LF", "M6"))

ggplot(Rich, aes(x=Filter_Type, y=Richness)) + geom_boxplot(aes(fill = Filter_Type), outlier.shape = NA) + geom_jitter(aes(shape = Date), col = "#b5b1b2", size = 3.5, width = 0.2) + theme_classic() + scale_fill_manual(values=c("#4363d8" ,"#000075","#808080","#cc0011"))


###Species accumulation curve (https://rdrr.io/rforge/vegan/man/specaccum.html). These were done seperately for bestperforming sampling setup and the rest of sampling setups.

#all other samplers
Tab3 <- read.csv2("Main_16S_Others.csv", row.names = 1, check.names = FALSE)
Tab4 <- Tab3[4:13]
sp1 <- specaccum(Tab4)
sp2 <- specaccum(Tab4, "random")
sp2
summary(sp2)
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", xlim=c(0, 20), ylim=c(0, 15))
boxplot(sp2, col="yellow", add=TRUE, pch="+")

##only for best sampler
Tab5 <- read.csv2("Main_16S_M6.csv", row.names = 1, check.names = FALSE)
Tab6 <- Tab5[4:15]
sp1 <- specaccum(Tab6)
sp2 <- specaccum(Tab6, "random")
sp2
summary(sp2)
plot(sp1, ci.type="poly", col="red", lwd=2, ci.lty=0, ci.col="lightgray", xlim=c(0, 20), ylim=c(0, 15))
boxplot(sp2, col="lightpink", add=TRUE, pch="+")


### Community level analyses. Here we are only using samples with taxa detections.
#Bring in the datasheet
Tab <- read.csv2("Main_16S.csv", row.names = 1, check.names = FALSE)

Tab2 <- Tab[4:17]
Meta <- Tab[1:3]

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

#download the coordinates and then add coordinates to the 16S_Stat document
write.csv2(princoor.pathotype.data, file = "PCOA_coord_16S_Filt_Raup.csv")

# Investigate the betadispersion (heterogeneity across the samples)
mod3 <- betadisper(dis,Meta$Filter_Type, type = "centroid", bias.adjust = TRUE)

#download the betadisperse output and add the data to the 16S_Stat document
write.csv2(mod3$distances, file = "Raup_beta_main_16S.csv")

#Conduct Permutational multivariate analyses tests (PERMANOVAS) using the adonis2 function. We use the same distance matrix that we used for PCOA.
adonis2(formula = dis~ Filter_Type, data = Meta, permutations = 10000, by = "margin", strata = Meta$Date)

#If PERMANOVA significant we can use the pairwise.adonis function to test pairwise differences
pairwise.adonis(dis, Meta$Filter_Type, p.adjust.m = "bonferroni", reduce = NULL,
                perm = 10000)

#Now bring in the stat document with the PCOA coordinates and betadispersal values
Data <- read.csv2("16S_Stat.csv", row.names = 1, check.names = FALSE)
Data$Filter_Type <- factor(Data$Filter_Type, levels = c("F8", "F7_FB", "F7_LF", "M6"))

#Plot the PCOA plot
ggplot(data = Data, aes(x = PCOA_X, y = PCOA_Y)) + geom_point(size=7, aes(colour = Filter_Type, shape = Date))+ scale_color_manual(values=c("#4363d8" ,"#000075","#808080", "#cc0011")) + theme_classic() + scale_shape_manual(values=c(16, 17, 15)) + stat_ellipse(aes(group = Filter_Type, col = Filter_Type), linetype = 2, type = "t", level = 0.80)

### box plots and linear mixed models of betadisperse
ggplot(Data, aes(x=Filter_Type, y=Beta_disp)) + geom_boxplot(aes(fill = Filter_Type), outlier.shape = NA) + geom_jitter(aes(shape = Date), col = "#b5b1b2", size = 3.5, width = 0.2) + theme_classic() + scale_fill_manual(values=c("#4363d8" ,"#000075","#808080","#cc0011"))

model_Beta <- lmer(Beta_disp ~  Filter_Type + (1|Date), Data)
summary(model_Beta)
anova(model_Beta)

Pair <- pairs(emmeans(model_Beta,  ~Filter_Type))
summary(Pair)







 







