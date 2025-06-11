###This r script is used to analyse the data from 16S mam primers for testing effect of sample storage method on accumulation of vertebrate DNA in the study "Advanced airborne eDNA sampling allows robust spatiotemporal characterisation of vertebrate communities"
## Experiment 2_(Æbelø)

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
library(UpSetR)

###Analysis of species richness. We are using all samples including samples without detections.
Rich <- read.csv2("16S_Meth_All.csv", row.names = 1, check.names = FALSE)

#Testing the effect of filter types on species richness using linear mixed models
model1 <- lmer(Richness ~  Storage_type + (1|Sampler) , Rich)

summary(model1)
anova(model1)

#Investigating the pairwise differences based on linear mixed models
Pair <- pairs(emmeans(model1,  ~Storage_type))
summary(Pair)

#Plotting the richness data with box plots
Rich$Storage_type <- factor(Rich$Storage_type, levels = c("Frozen", "LB", "Buffer"))

#Geom points with merging pairs
ggplot(Rich, aes(x=Storage_type, y=Richness)) + geom_boxplot(aes(fill = Storage_type), outlier.shape = NA)+ geom_point(aes(size = 7))+ scale_fill_manual(values=c("#91de45" ,"#bc4ec2","#63d6e6")) + scale_color_manual(c("#000")) +theme_classic()+ geom_line(aes(group=Sampler), linetype = 2)

###Generating Upset plots
Up <- read.csv2("16S_Upset.csv", row.names = 1, check.names = FALSE)
Up <- as.data.frame(Up)
upset(Up, sets = c("Buffer","LB","Frozen" ), sets.bar.color = "#56B4E9",
      order.by = "degree", empty.intersections = "off", keep.order = TRUE)


###Species accumulation curve (https://rdrr.io/rforge/vegan/man/specaccum.html). These were done seperately for bestperforming sampling setup and the rest of sampling setups.
#All other samples
Tab5 <- read.csv2("Storage_16S_Others.csv", row.names = 1, check.names = FALSE)
Tab6 <- Tab5[3:9]
sp1 <- specaccum(Tab6)
sp2 <- specaccum(Tab6, "random")
sp2
summary(sp2)
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",xlim=c(0, 35), ylim=c(0, 15))
boxplot(sp2, col="yellow", add=TRUE, pch="+")

##only for best sampler
Tab3 <- read.csv2("Storage_16S_dry.csv", row.names = 1, check.names = FALSE)
Tab4 <- Tab3[3:14]
sp1 <- specaccum(Tab4)
sp2 <- specaccum(Tab4, "random")
sp2
summary(sp2)
plot(sp1, ci.type="poly", col="red", lwd=2, ci.lty=0, ci.col="lightgray", xlim=c(0, 35), ylim=c(0, 15))
boxplot(sp2, col="lightpink", add=TRUE, pch="+")

### Community level analyses. Here we are only using samples with taxa detections.
#Bring in the datasheet
Tab <- read.csv2("Storage_16S.csv", row.names = 1, check.names = FALSE)

Tab2 <- Tab[3:16]
Meta <- Tab[1:2]

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

#download the coordinates and then add coordinates to the 16S_Meth_Stat document
write.csv2(princoor.pathotype.data, file = "16S_Method_Cor_Raup.csv")

#Conduct Permutational multivariate analyses tests (PERMANOVAS) using the adonis2 function. We use the same distance matrix that we used for PCOA.
adonis2(formula = dis~ Storage_type, data = Meta, permutations = 10000, by = "margin", strata = Meta$Sampler)

#If PERMANOVA significant we can use the pairwise.adonis function to test pairwise differences
pairwise.adonis(dis, Meta$Storage_type, p.adjust.m = "bonferroni", reduce = NULL,
                perm = 10000)


#Now bring in the stat document with the PCOA coordinates
Data <- read.csv2("16S_Meth_Stat.csv", row.names = 1, check.names = FALSE)
Data$Storage_type <- factor(Data$Storage_type, levels = c("Frozen", "LB", "Buffer"))

#Plot the PCOA plot
ggplot(data = Data, aes(x = Cor_X, y = Cor_Y)) + geom_point(size=7, aes(colour = Storage_type))+ scale_color_manual(values=c("#91de45" ,"#bc4ec2","#63d6e6")) + theme_classic() + scale_shape_manual(values=c(16, 17, 15)) + geom_line(aes(group=Sampler))+ stat_ellipse(aes(group = Storage_type, col = Storage_type), linetype = 2,lwd = 1.2, type = "t", level = 0.95)
ggplot(data = Data, aes(x = Cor_X, y = Cor_Y)) + geom_point(size=7, aes(colour = Storage_type))+ scale_color_manual(values=c("#91de45" ,"#bc4ec2","#63d6e6")) + theme_classic() + scale_shape_manual(values=c(16, 17, 15)) + geom_text(aes(label = Sampler))


#####Mantel tests###
#make sure sampler names are in matching order between OTU tables and the coordiante tables
#Getting distance matrix between samplers
library(sf)

#Bring in the coordinate
Frozen <-read.csv2("Coordinates_Frozen.csv", row.names = 1, check.names = FALSE)
LB<-read.csv2("Coordinates_LB.csv", row.names = 1, check.names = FALSE)

species_sf <- st_as_sf(Frozen, coords = c("Longitude", "Latitude"), crs = 4326)
Frozen_m <- st_distance(species_sf)
LB_m <- st_distance(species_sf)

write.csv2(as.matrix(Frozen_m), file = "Frozen_distance_m.csv")
write.csv2(as.matrix(LB_m), file = "LB_distance_m.csv")

# Bring in the OTU tables
Frozen_OTU <-read.csv2("Frozen_Mantel_OTU.csv", row.names = 1, check.names = FALSE)
LB_OTU<-read.csv2("LB_Mantel_OTU.csv", row.names = 1, check.names = FALSE)

#Calculate community distances between samples
Frozen_com_dist <- raupcrick(Frozen_OTU, null = "r1", nsimul = 999)
LB_com_dist <- raupcrick(LB_OTU, null = "r1", nsimul = 999)

write.csv2(as.matrix(Frozen_com_dist), file = "Frozen_OTU_dist_Raup.csv")
write.csv2(as.matrix(LB_com_dist), file = "LB_OTU_dist_Raup.csv")

#running the mantel test
Man<- mantel(LB_com_dist, LB_m, method="pearson", permutations=10000, strata = NULL,na.rm = FALSE)
Man

#Now plot the two distance matrixes together
Dist <-read.csv2("Mantel_figure.csv", row.names = 1, check.names = FALSE)

ggplot(Dist, aes(x=Dist_m, y=OTU_dist, col = Storage_type))+  geom_point(aes(col = Storage_type, size = 7))+ scale_color_manual(values=c("#91de45" ,"#bc4ec2"))+theme_classic()+ geom_smooth(method=lm , se=TRUE)



