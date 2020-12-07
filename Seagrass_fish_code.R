
library(car)
library(lme4)
library(VIF)
library(ggplot2)
library(ggfortify)
library(dplyr)
library(GGally)
library(CCA)
require(MuMIn)
library(vegan)
library(PerformanceAnalytics)
library(rcompanion)
require(nlme)
library(visreg)
library(EnvStats)
library(viridis)
library(wesanderson)
######ALL SEAGRASS VARIABLES#####

# Set working directory (Dropbox R Scripts)
setwd("~/Dropbox/Phd/Paper II - Seagrass diversity fisheries/Data")

# load seagrass variables into R and summarise (all data)
seagrass_var<-read.csv("seagrass_data.csv", header = TRUE, sep=";")
as.matrix(seagrass_var)
head(seagrass_var)
names(seagrass_var)
summary(seagrass_var)
seagrass_var$site<-as.factor(seagrass_var$site)

scatterplot(Zanz_fish2$land.use.pc1 ~ seagrass_var_average$X._epiphyte)


# compute PCA for all seagrass traits and decompose to creat coefficient of variation
traits_all<-prcomp(~shoot_density+canopy+leaf+width+no_leaves1, data = seagrass_var, center=TRUE, scale=TRUE)
autoplot(prcomp(~shoot_density+canopy+leaf+width+no_leaves1, data = seagrass_var, center=TRUE, scale = TRUE), data = seagrass_var, colour = 'site',
         loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 4, loadings.label.colour = "black") + scale_color_viridis(discrete = TRUE, option = "D")

summary(traits_all)
traits_ev <- traits_all$sdev^2
traits_ev
prcomp(~shoot_density+canopy+leaf+width+no_leaves1, data = seagrass_var, center=TRUE, scale=TRUE)

traits.data_all<-cbind(seagrass_var[c(1:2)], predict(traits_all, newdata = seagrass_var))

require(plyr)
traits.data_all<-ddply(traits.data_all, .(sample), summarize,  trait.mean.pc1=mean(PC1), trait.sd.pc1=sd(PC1))

traits.data_all<-ddply(traits.data_all, .(sample), summarize,  trait.co.var.pc1=(trait.sd.pc1/trait.mean.pc1)*100)

#Average all seagrass variables
seagrass_var_average<- seagrass_var %>%
  group_by(sample,site) %>% 
  summarise_all(mean)

seagrass_var_average<-data.frame(seagrass_var_average)
str(seagrass_var_average)

skew<-skewness(seagrass_var_average[c(3:23)])

# compute PCA for composition/cover variables and decompose
composition<-prcomp(~Cr+Cs+Hu+Ho+Th+Tc+Ea+Si+Hs, data = seagrass_var_average, center=TRUE, scale=FALSE)
autoplot(prcomp(~Cr+Cs+Hu+Ho+Th+Tc+Ea+Si+Hs, data = seagrass_var_average, center=TRUE, scale = FALSE), data = seagrass_var_average, colour = 'site',
         loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 4, loadings.label.colour = "black") + scale_color_viridis(discrete = TRUE, option = "D")
summary(composition)
composition_ev <- composition$sdev^2
composition_ev
prcomp(~Cr+Cs+Hu+Ho+Th+Tc+Ea+Si+Hs, data = seagrass_var_average, center=TRUE, scale=FALSE)

composition.data <- cbind(seagrass_var_average[c(1:2)], predict(composition, newdata = seagrass_var_average))
head(composition.data, 5)
require(plyr)
composition.data<-ddply(composition.data, .(sample), summarize,  comp.pc1=mean(PC1), comp.pc2=mean(PC2))

write.table(composition.data, "~/Dropbox/Phd/Paper II - Seagrass diversity fisheries/Data/pca/composition.data.txt", sep="\t")

# compute PCA for seagrass traits and decompose
traits<-prcomp(~shoot_density+canopy+leaf+width+no_leaves1, data = seagrass_var_average, center=TRUE, scale=TRUE)
autoplot(prcomp(~shoot_density+canopy+leaf+width+no_leaves1, data = seagrass_var_average, center=TRUE, scale = TRUE), data = seagrass_var_average, colour = 'site',
         loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 4, loadings.label.colour = "black") + scale_color_viridis(discrete = TRUE, option = "D")
summary(traits)
traits_ev <- traits$sdev^2
traits_ev
prcomp(~shoot_density+canopy+leaf+width+no_leaves1, data = seagrass_var_average, center=TRUE, scale=TRUE)

traits.data<-cbind(seagrass_var_average[c(1:2)], predict(traits, newdata = seagrass_var_average))
require(plyr)
traits.data<-ddply (traits.data, .(sample), summarize,  trait.pc1=mean(PC1), trait.pc2=mean(PC2))

write.table(traits.data, "~/Dropbox/Phd/Paper II - Seagrass diversity fisheries/Data/traits.data.txt", sep="\t")


#####MAIN VARIABLES AND PREDICTORS #######
# load seagrass and fish variables into R and summarise (avergae)
setwd("~/Dropbox/Phd/Paper II - Seagrass diversity fisheries/Data")
Zanz_fish<-read.csv("enviro_seagrass_fish2.csv", header = TRUE, sep=";")
as.matrix(Zanz_fish)
names(Zanz_fish)
head(Zanz_fish)
summary(Zanz_fish)
Zanz_fish$Fishing.pressure<-factor(Zanz_fish$Fishing.pressure, c("low", "high"), levels = c("low", "high"))
Zanz_fish$Fishing.pressure<-as.numeric(Zanz_fish$Fishing.pressure)
skew_bigdata<-PerformanceAnalytics:: skewness(Zanz_fish[c(3:24)])

MaxN_average<- Zanz_fish %>%
  group_by(Site) %>% 
  summarise_all(mean)


#transform land use  Bareground+Human.development+Aquaculture+Farmland+Vegetation
Zanz_fish$Bareground<-sqrt(Zanz_fish$Bareground)
Zanz_fish$Human.development<-sqrt(Zanz_fish$Human.development)
Zanz_fish$Aquaculture<-sqrt(Zanz_fish$Aquaculture)
Zanz_fish$Farmland<-sqrt(Zanz_fish$Farmland)
Zanz_fish$Vegetation<-sqrt(Zanz_fish$Vegetation)
Zanz_fish$Ocean<-sqrt(Zanz_fish$Ocean)

# compute PCA for landuse and decompose + bind into dataframe
land.use<-prcomp (~Bareground+Human.development+Aquaculture+Farmland+Vegetation+D.catchment+Ocean, data = Zanz_fish, center=TRUE, scale=TRUE)
autoplot(prcomp(~Bareground+Human.development+Aquaculture+Farmland+Vegetation+D.catchment+Ocean, data = Zanz_fish, center=TRUE, scale =TRUE), data = Zanz_fish, colour = 'Site',
         loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 4, loadings.label.colour = "black")+scale_color_viridis(discrete = TRUE, option = "D")
summary(land.use)
prcomp(~Bareground+Human.development+Aquaculture+Farmland+Vegetation+D.catchment+Ocean, data = Zanz_fish, center=TRUE, scale=TRUE)

# Bind new dataset including land use, composition data, traits and traits co.var
Zanz_fish<- cbind(Zanz_fish,land.use.pc1=predict(land.use)[,1],composition.data[c(2)],traits.data[c(2)], traits.data_all[c(2)])

# invert land use and traits
Zanz_fish$land.use.pc1<-Zanz_fish$land.use.pc1 * -1
Zanz_fish$trait.pc1<-Zanz_fish$trait.pc1 * -1


# Make new dataset removing old variables
Zanz_fish1<-data.frame(Zanz_fish[c(2:13,21:28)])

chart.Correlation(Zanz_fish1[c(2:20)], histogram=TRUE, pch="+")
skew_bigdata<-PerformanceAnalytics:: skewness(Zanz_fish1[c(2:20)])

library(ggpubr)

fishingplot<-ggplot(Zanz_fish1, aes(x=Fishing.pressure, y=land.use.pc1)) + 
  geom_point(color="#1B9E77")+
  geom_smooth(method=lm, color="#1B9E77", fill="#1B9E77")+
  xlab("Fishing Pressure")+
  ylab("Land-use")+
  theme_gray()+
  stat_cor(method = "pearson", label.x.npc=0.5, p.accuracy = 0.001)

mangroveplot<-ggplot(Zanz_fish1, aes(x=D.mangrove, y=land.use.pc1)) + 
  geom_point(color="#1B9E77")+
  geom_smooth(method=lm, color="#1B9E77", fill="#1B9E77")+
  xlab("Distance to Mangrove")+
  ylab("Land-use")+
  theme_gray()+
  stat_cor(method = "pearson", label.x.npc=0.5, p.accuracy = 0.001)

coralplot<-ggplot(Zanz_fish1, aes(x=D.coral, y=land.use.pc1)) + 
  geom_point(color="#1B9E77")+
  geom_smooth(method=lm, color="#1B9E77", fill="#1B9E77")+
  xlab("Distance to Coral")+
  ylab("Land-use")+
  theme()+
  stat_cor(method = "pearson", label.y=5, p.accuracy = 0.001)

coralplot
Landingplot<-ggplot(Zanz_fish1, aes(x=D.landing.site, y=land.use.pc1)) + 
  geom_point(color="#1B9E77")+
  geom_smooth(method=lm, color="#1B9E77", fill="#1B9E77")+
  xlab("Distance to Landing Site")+
  ylab("Land-use")+
  theme_gray()+
  stat_cor(method = "pearson", label.x.npc=0.5, p.accuracy = 0.001)

LandUseplot <- cowplot::plot_grid(fishingplot, mangroveplot, coralplot, Landingplot, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), align = "hv")

LandUseplot



# transform skewed data
Zanz_fish1$MaxN<- log1p(Zanz_fish1$MaxN)
Zanz_fish1$Sric <- log(Zanz_fish1$Sric)
Zanz_fish1$FRic <- log(Zanz_fish1$FRic)
Zanz_fish1$FDis <- log1p(Zanz_fish1$FDis)
Zanz_fish1$Depth <- log(Zanz_fish1$Depth)
Zanz_fish1$D.coral <- log(Zanz_fish1$D.coral)
Zanz_fish1$TM.MaxN<- log1p(Zanz_fish1$TM.MaxN)
Zanz_fish1$SST.MaxN<- log1p(Zanz_fish1$SST.MaxN)
Zanz_fish1$LH.MaxN<-log1p(Zanz_fish1$LH.MaxN)
Zanz_fish1$LV.MaxN<- log1p(Zanz_fish1$LV.MaxN)
Zanz_fish1$comp.pc1 <- sign(Zanz_fish1$comp.pc1) * (abs(Zanz_fish1$comp.pc1))^(1/3)
Zanz_fish1$trait.co.var.pc1 <- sign(Zanz_fish1$trait.co.var.pc1) * (abs(Zanz_fish1$trait.co.var.pc1))^(1/3)

# Make new dataset removing old variables and scale
Zanz_fish2<-data.frame(Zanz_fish1[c(1:20)])
                  
# correlation matrix for new variarbles
chart.Correlation(Zanz_fish2[c(2:20)], histogram=TRUE, pch="+") #D.mangrove, D.coral, D.landing site, Fishing pressure correlated with Land Use.
skew_bigdata<-PerformanceAnalytics:: skewness(Zanz_fish2[c(2:20)])
names(Zanz_fish2)

ggplot(Zanz_fish2, aes(x=trait.pc1, y=Depth)) + 
  geom_point(color="#1B9E77")+
  geom_smooth(method=lm, color="#1B9E77", fill="#1B9E77")+
  theme_gray()


#Build model using all predictors and to detect multicollinearity
mod1 <- lme(MaxN~Sric+FDis+FRic+SCov+Depth+land.use.pc1+
              comp.pc1+trait.pc1+trait.co.var.pc1, 
             random=~1|Site, data = Zanz_fish2, method = "REML")

mod2 <- lme(Fish.Ric~Sric+FDis+FRic+SCov+Depth+land.use.pc1+
              comp.pc1+trait.pc1+trait.co.var.pc1, 
            random=~1|Site, data = Zanz_fish2, method = "REML")

car::vif(mod1) #all variables have good VIF
car::vif(mod2) #all variables have good VIF

###### Hypothesis models MaxN #####
#seagrass cover
require(lme4)
require(effects)
require(ggplot2)
Seagrass_Cover <- lmer(MaxN~SCov*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Cover, test="Chisq")
summary(Seagrass_Cover, correlation=FALSE)


#Seagrass Richness
Seagrass_Richness <- lmer(MaxN~Sric*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1 (Seagrass_Richness, test="Chisq") #remove interaction
Seagrass_Richness <- lmer(MaxN~Sric+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Richness, test="Chisq") #remove land-use
Seagrass_Richness <- lmer(MaxN~Sric+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Richness, test="Chisq") #could not be improved without removing Seagrass Richness
summary(Seagrass_Richness, correlation=FALSE)

#Functional dispersion
Functional_Dispersion <- lmer(MaxN~FDis*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion, test="Chisq") #remove interaction
Functional_Dispersion <- lmer(MaxN~FDis+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion, test="Chisq") 
summary(Functional_Dispersion, correlation=FALSE)

#Functional Richness
Functional_Richness <- lmer(MaxN~FRic*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness, test="Chisq") #remove interaction
Functional_Richness <- lmer(MaxN~FRic+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness, test="Chisq") #remove land-use
Functional_Richness <- lmer(MaxN~FRic+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness, test="Chisq") #could not be improved without removing Functional Richness
summary(Functional_Richness, correlation=FALSE)

#Meadow Trtaits
Meadow_Traits <- lmer(MaxN~trait.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits, test="Chisq") #drop interaction
Meadow_Traits <- lmer(MaxN~trait.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits, test="Chisq") #drop land use
Meadow_Traits <- lmer(MaxN~trait.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits, test="Chisq")
summary(Meadow_Traits, correlation=FALSE)

library(nlme)
Meadow_Traits <- lme(MaxN~trait.pc1+Depth, random = ~ 1 | Site, data=Zanz_fish2, method = "ML")

require(visreg)
visreg(Meadow_Traits, "trait.pc1", band=TRUE)

require(jtools)
Traits_plot<-effect_plot(Meadow_Traits, pred = "trait.pc1", interval = TRUE, plot.points = TRUE, colors = "#1B9E77", point.color = "#1B9E77") +xlab("Meadow Stucture") +
  ylab("Log Fish Abundance")+theme(legend.position = "none", text = element_text(size=14))+coord_cartesian(ylim = c(-1, 6))
Traits_plot

Depth_plot <-effect_plot(Meadow_Traits, pred = "Depth", interval = TRUE, plot.points = TRUE, colors = "#1B9E77", point.color = "#1B9E77") +xlab("Depth") +
  ylab("")+theme(legend.position = "none", text = element_text(size=14))+coord_cartesian(ylim = c(-1, 6))

library(egg)
MaxNPlotModel<-ggarrange(Traits_plot, Depth_plot, ncol = 2, nrow = 1, labels = c("A", "B"))

MaxNPlotModel <- cowplot::plot_grid(Traits_plot, Depth_plot, ncol=2, labels = c("A", "B"), align = "hv")
MaxNPlotModel

#Seagrass Composition
Seagrass_Composition <- lmer(MaxN~comp.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition, test="Chisq") #Drop interaction
Seagrass_Composition <- lmer(MaxN~comp.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition, test="Chisq") #drop land-use
Seagrass_Composition <- lmer(MaxN~comp.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition, test="Chisq") 
summary(Seagrass_Composition, correlation=FALSE)

#Meadow trait variability
Meadow_Trait_Variability <- lmer(MaxN~trait.co.var.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability) #drop interaction
Meadow_Trait_Variability <- lmer(MaxN~trait.co.var.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability) #drop land-use
Meadow_Trait_Variability <- lmer(MaxN~trait.co.var.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability)
summary(Meadow_Trait_Variability, correlation=FALSE)

#null model
null_model <- lmer(MaxN~land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
summary(null_model, correlation=FALSE)

require(AICcmodavg)
MaxN_model_list<-list(null_model, Seagrass_Richness, Seagrass_Cover, Functional_Dispersion, Functional_Richness, Meadow_Traits, Seagrass_Composition, Meadow_Trait_Variability)
model_names<-c("Null Model", "Seagrass Richness", "Seagrass Cover", "Functional Dispersion", "Functional Richness", "Meadow Traits", "Seagrass Composition", "Meadow Trait Variability")
modelsel<-aictab(MaxN_model_list, model_names, second.ord=T)
modelsel
write.table(modelsel, "~/Dropbox/Phd/Paper II - Seagrass diversity fisheries/Data/MaxN_Model_Selection.txt", sep="\t")

require(piecewiseSEM)
rsquared(MaxN_model_list)

anova(null_model, Seagrass_Richness, Seagrass_Cover, Functional_Dispersion, Functional_Richness, Meadow_Traits, Seagrass_Composition, Meadow_Trait_Variability)




###### Hypothesis models Richness #####
#seagrass cover
require(lme4)
require(effects)
require(ggplot2)
Seagrass_Cover_Ric <- lmer(Fish.Ric~SCov*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Cover_Ric, test="Chisq")
summary(Seagrass_Cover_Ric, correlation=FALSE)


visreg(Seagrass_Cover_Ric, "SCov", by="land.use.pc1", overlay = TRUE, gg=TRUE, ylab="Fish Richness", xlab="Seagrass Cover")

Cover_Plot <- visreg(Seagrass_Cover_Ric, "SCov", by="land.use.pc1", overlay = TRUE, gg=TRUE, ylab="Fish Richness", xlab="Seagrass Cover") + ylim(-5, 15) + theme(text = element_text(size = 14)) + theme(legend.position=c(0.5,0.9))
Depth_plot_richness <- visreg(Seagrass_Cover_Ric, "Depth", gg=TRUE, ylab="", xlab="Log Depth") + ylim(-5, 15) + theme(text = element_text(size = 14))

RichnessPlotModel<-egg:: ggarrange(Cover_Plot, Depth_plot_richness, ncol = 2, nrow = 1)

library(interactions)
Cover_Plot<- interact_plot(Seagrass_Cover_Ric, pred = "SCov", modx = land.use.pc1, interval = TRUE, partial.residuals = TRUE, colors = "Dark2",) +xlab("Seagrass Cover") +
  ylab("Fish Richness")+theme(text = element_text(size=14), legend.position = c(0.7, 0.2),
                              legend.direction = "vertical")+coord_cartesian(ylim = c(-5, 15))
Cover_Plot


get_colors("Dark2")

Depth_plot_richness <- effect_plot(Seagrass_Cover_Ric, interval=TRUE, pred = Depth, partial.residuals = TRUE, colors = "#1B9E77", point.color = "#1B9E77") +xlab("Depth") +
  ylab("")+theme(text = element_text(size=14))+coord_cartesian(ylim = c(-5, 15))
Depth_plot_richness

RichnessPlotModel<-egg:: ggarrange(Cover_Plot, Depth_plot_richness, ncol = 2, nrow = 1)
RichnessPlotModel <- cowplot::plot_grid(Cover_Plot, Depth_plot_richness, ncol=2, nrow = 1, labels = c("A", "B"), align = "hv")
RichnessPlotModel

coord_cartesian(ylim = c(0, 10))
#Seagrass Richness
Seagrass_Richness_Ric <- lmer(Fish.Ric~Sric*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Richness_Ric, test="Chisq") #remove interaction
Seagrass_Richness_Ric <- lmer(Fish.Ric~Sric+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Richness_Ric, test="Chisq") #could not be improved without removing Seagrass Richness
summary(Seagrass_Richness_Ric, correlation=FALSE)

#Functional dispersion
Functional_Dispersion_Ric <- lmer(Fish.Ric~FDis*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion_Ric, test="Chisq") #drop interaction
Functional_Dispersion_Ric <- lmer(Fish.Ric~FDis+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion_Ric, test="Chisq")
summary(Functional_Dispersion_Ric, correlation=FALSE)

#Functional Richness
Functional_Richness_Ric <- lmer(Fish.Ric~FRic*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness_Ric, test="Chisq") #remove interaction
Functional_Richness_Ric <- lmer(Fish.Ric~FRic+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness_Ric, test="Chisq") #could not be improved without removing Functional Richness
summary(Functional_Richness_Ric, correlation=FALSE)

#Meadow Trtaits
Meadow_Traits_Ric <- lmer(Fish.Ric~trait.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_Ric, test="Chisq") #drop interaction
Meadow_Traits_Ric <- lmer(Fish.Ric~trait.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_Ric, test="Chisq") #drop land use
Meadow_Traits_Ric <- lmer(Fish.Ric~trait.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_Ric, test="Chisq")
summary(Meadow_Traits_Ric, correlation=FALSE)

visreg(Meadow_Traits_Ric, "trait.pc1", by="land.use.pc1", overlay = TRUE, gg=TRUE, ylab="Fish Richness", xlab="Meadow Structure")


Trait_Plot_Rich <- visreg(Meadow_Traits_Ric, "trait.pc1", gg=TRUE, ylab="Fish Richness", xlab="Meadow Structure") + ylim(-5, 15) + theme(text = element_text(size = 14))
Depth_plot_Rich <- visreg(Meadow_Traits_Ric, "Depth", gg=TRUE, ylab="", xlab="Log Depth") + ylim(-5, 15) + theme(text = element_text(size = 14))

RichnessTraitPlotModel<-egg:: ggarrange(Trait_Plot_Rich, Depth_plot_Rich, ncol = 2, nrow = 1)

require(jtools)
Trait_Plot_Rich<-effect_plot(Meadow_Traits_Ric, pred = "trait.pc1", interval = TRUE, plot.points = TRUE, colors = "#1B9E77", point.color = "#1B9E77") +xlab("Meadow Stucture") +
  ylab("Fish Richness")+theme(legend.position = "none", text = element_text(size=14))+coord_cartesian(ylim = c(-5, 15))
Trait_Plot_Rich

Depth_plot_Rich <-effect_plot(Meadow_Traits_Ric, pred = "Depth", interval = TRUE, plot.points = TRUE, colors = "#1B9E77", point.color = "#1B9E77") +xlab("Depth") +
  ylab("")+theme(legend.position = "none", text = element_text(size=14))+coord_cartesian(ylim = c(-5, 15))

library(egg)
RichnessPlotModel<-ggarrange(Trait_plot_Rich, Depth_plot_Rich, ncol = 2, nrow = 1, labels = c("A", "B"))

RichnessPlotModel <- cowplot::plot_grid(Trait_Plot_Rich, Depth_plot_Rich, ncol=2, labels = c("A", "B"), align = "hv")
RichnessPlotModel


#Seagrass Composition
Seagrass_Composition_Ric <- lmer(Fish.Ric~comp.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition_Ric, test="Chisq") #Drop interaction
Seagrass_Composition_Ric <- lmer(Fish.Ric~comp.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition_Ric, test="Chisq") #could not be improved without removing Seagrass composition
summary(Seagrass_Composition_Ric, correlation=FALSE)

#Meadow trait variability
Meadow_Trait_Variability_Ric <- lmer(Fish.Ric~trait.co.var.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability_Ric, test="Chisq") #drop interaction
Meadow_Trait_Variability_Ric <- lmer(Fish.Ric~trait.co.var.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability_Ric, test="Chisq")
summary(Meadow_Trait_Variability_Ric, correlation=FALSE)

#null model
null_model_Ric <- lmer(Fish.Ric~land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
summary(null_model_Ric, correlation=FALSE)

#compare models
require(AICcmodavg)
Ric_model_list<-list(null_model_Ric, Seagrass_Richness_Ric, Seagrass_Cover_Ric, Functional_Dispersion_Ric, Functional_Richness_Ric, Meadow_Traits_Ric, Seagrass_Composition_Ric, Meadow_Trait_Variability_Ric)
model_names<-c("Null Model", "Seagrass Richness", "Seagrass Cover", "Functional Dispersion", "Functional Richness", "Meadow Traits", "Seagrass Composition", "Meadow Trait Variability")
modelsel<-aictab(Ric_model_list, model_names, second.ord=T)
write.table(modelsel, "~/Dropbox/Phd/Paper II - Seagrass diversity fisheries/Data/Ric_Model_Selection.txt", sep="\t")
modelsel
require(piecewiseSEM)
rsquared(Ric_model_list)

anova(null_model_Ric, Seagrass_Richness_Ric, Seagrass_Cover_Ric, Functional_Dispersion_Ric, Functional_Richness_Ric, Meadow_Traits_Ric, Seagrass_Composition_Ric, Meadow_Trait_Variability_Ric)

#######Hypothesis Fish Value######
#Seagrass Cover
require(lme4)
require(effects)
require(ggplot2)
Seagrass_Cover_TM <- lmer(TM.MaxN~SCov*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Cover_TM, test="Chisq") #Drop interaction
Seagrass_Cover_TM <- lmer(TM.MaxN~SCov+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Cover_TM, test="Chisq") 
summary(Seagrass_Cover_TM, correlation=FALSE)

Seagrass_Cover_SST <- lmer(SST.MaxN~SCov*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Cover_SST, test="Chisq") #Drop interaction
Seagrass_Cover_SST <- lmer(SST.MaxN~SCov+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Cover_SST, test="Chisq") 
summary(Seagrass_Cover_SST, correlation=FALSE)

Seagrass_Cover_LH <- lmer(LH.MaxN~SCov*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Cover_LH, test="Chisq") #Drop interaction
Seagrass_Cover_LH <- lmer(LH.MaxN~SCov+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Cover_LH, test="Chisq")
summary(Seagrass_Cover_LH, correlation=FALSE)

Seagrass_Cover_LV <- lmer(LV.MaxN~SCov*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Cover_LV, test="Chisq") 
summary(Seagrass_Cover_LV, correlation=FALSE)

#Seagrass Richness

Seagrass_Richness_TM <- lmer(TM.MaxN~Sric*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Richness_TM, test="Chisq") #Drop interaction
Seagrass_Richness_TM <- lmer(TM.MaxN~Sric+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Richness_TM, test="Chisq")
summary(Seagrass_Richness_TM, correlation=FALSE)

Seagrass_Richness_SST <- lmer(SST.MaxN~Sric*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Richness_SST, test="Chisq") #Drop interaction
Seagrass_Richness_SST <- lmer(SST.MaxN~Sric+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Richness_SST, test="Chisq") #drop land use
Seagrass_Richness_SST <- lmer(SST.MaxN~Sric+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Richness_SST, test="Chisq")
summary(Seagrass_Richness_SST, correlation=FALSE)

Seagrass_Richness_LH <- lmer(LH.MaxN~Sric*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Richness_LH, test="Chisq") #Drop interaction
Seagrass_Richness_LH <- lmer(LH.MaxN~Sric+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Richness_LH, test="Chisq") #no improvment without loss of seagrass variable
summary(Seagrass_Richness_LH, correlation=FALSE)

Seagrass_Richness_LV <- lmer(LV.MaxN~Sric*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Richness_LV, test="Chisq")
summary(Seagrass_Richness_LV, correlation=FALSE)

#Functional Dispersion

Functional_Dispersion_TM <- lmer(TM.MaxN~FDis*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion_TM, test="Chisq") #Drop interaction
Functional_Dispersion_TM <- lmer(TM.MaxN~FDis+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion_TM, test="Chisq")
summary(Functional_Dispersion_TM, correlation=FALSE)

Functional_Dispersion_SST <- lmer(SST.MaxN~FDis*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion_SST, test="Chisq") #Drop interaction
Functional_Dispersion_SST <- lmer(SST.MaxN~FDis+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion_SST, test="Chisq")
summary(Functional_Dispersion_SST, correlation=FALSE)

Functional_Dispersion_LH <- lmer(LH.MaxN~FDis*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion_LH, test="Chisq") #Drop interaction
Functional_Dispersion_LH <- lmer(LH.MaxN~FDis+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion_LH, test="Chisq")
summary(Functional_Dispersion_LH, correlation=FALSE)

visreg(Functional_Dispersion_LH)

Functional_Dispersion_LV <- lmer(LV.MaxN~FDis*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion_LV, test="Chisq") #Drop interaction
Functional_Dispersion_LV <- lmer(LV.MaxN~FDis+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion_LV, test="Chisq") #Drop land-use
Functional_Dispersion_LV <- lmer(LV.MaxN~FDis+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Dispersion_LV, test="Chisq")
summary(Functional_Dispersion_LV, correlation=FALSE)

#Functional Richness

Functional_Richness_TM <- lmer(TM.MaxN~FRic*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness_TM, test="Chisq") #Drop interaction
Functional_Richness_TM <- lmer(TM.MaxN~FRic+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness_TM, test="Chisq") #no improvment without loss of seagrass variable
summary(Functional_Richness_TM, correlation=FALSE)

Functional_Richness_SST <- lmer(SST.MaxN~FRic*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness_SST, test="Chisq") #Drop interaction
Functional_Richness_SST <- lmer(SST.MaxN~FRic+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness_SST, test="Chisq") #drop land use
Functional_Richness_SST <- lmer(SST.MaxN~FRic+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness_SST, test="Chisq")
summary(Functional_Richness_SST, correlation=FALSE)

Functional_Richness_LH <- lmer(LH.MaxN~FRic*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness_LH, test="Chisq") #drop interaction
Functional_Richness_LH <- lmer(LH.MaxN~FRic+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness_LH, test="Chisq")
summary(Functional_Richness_LH, correlation=FALSE)

Functional_Richness_LV <- lmer(LV.MaxN~FRic*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness_LV, test="Chisq") #drop interaction
Functional_Richness_LV <- lmer(LV.MaxN~FRic+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness_LV, test="Chisq") #drop land use
Functional_Richness_LV <- lmer(LV.MaxN~FRic+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Functional_Richness_LV, test="Chisq")
summary(Functional_Richness_LV, correlation=FALSE)


#Meadow Traits

Meadow_Traits_TM <- lmer(TM.MaxN~trait.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_TM, test="Chisq") #Drop interaction
Meadow_Traits_TM <- lmer(TM.MaxN~trait.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_TM, test="Chisq") #drop land-use
Meadow_Traits_TM <- lmer(TM.MaxN~trait.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_TM, test="Chisq")
summary(Meadow_Traits_TM, correlation=FALSE)

Meadow_Traits_SST <- lmer(SST.MaxN~trait.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_SST, test="Chisq") #Drop interaction
Meadow_Traits_SST <- lmer(SST.MaxN~trait.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_SST, test="Chisq") #drop land use
Meadow_Traits_SST <- lmer(SST.MaxN~trait.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_SST, test="Chisq")
summary(Meadow_Traits_SST, correlation=FALSE)

Meadow_Traits_LH <- lmer(LH.MaxN~trait.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_LH, test="Chisq") #Drop interaction
Meadow_Traits_LH <- lmer(LH.MaxN~trait.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_LH, test="Chisq") 
summary(Meadow_Traits_LH, correlation=FALSE)

Meadow_Traits_LV <- lmer(LV.MaxN~trait.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_LV, test="Chisq") #drop interaction
Meadow_Traits_LV <- lmer(LV.MaxN~trait.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_LV, test="Chisq") #drop land use
Meadow_Traits_LV <- lmer(LV.MaxN~trait.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Traits_LV, test="Chisq") 
summary(Meadow_Traits_LH, correlation=FALSE)

#Seagrass Composition
Seagrass_Composition_TM <- lmer(TM.MaxN~comp.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition_TM, test="Chisq") #drop interaction
Seagrass_Composition_TM <- lmer(TM.MaxN~comp.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition_TM, test="Chisq") #drop land use
Seagrass_Composition_TM <- lmer(TM.MaxN~comp.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition_TM, test="Chisq")
summary(Seagrass_Composition_TM, correlation=FALSE)

Seagrass_Composition_SST <- lmer(SST.MaxN~comp.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition_SST, test="Chisq") #drop interaction
Seagrass_Composition_SST <- lmer(SST.MaxN~comp.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition_SST, test="Chisq") #drop land use
Seagrass_Composition_SST <- lmer(SST.MaxN~comp.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition_SST, test="Chisq") 
summary(Seagrass_Composition_SST, correlation=FALSE)

Seagrass_Composition_LH <- lmer(LH.MaxN~comp.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition_LH, test="Chisq") #Drop interaction
Seagrass_Composition_LH <- lmer(LH.MaxN~comp.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition_LH, test="Chisq") #no improvment without loss of seagrass variable
summary(Seagrass_Composition_LH, correlation=FALSE)

Seagrass_Composition_LV <- lmer(LV.MaxN~comp.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Seagrass_Composition_LV, test="Chisq")
summary(Seagrass_Composition_LV, correlation=FALSE)

#Meadow Trait Variability

Meadow_Trait_Variability_TM <- lmer(TM.MaxN~trait.co.var.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability_TM, test="Chisq") #Drop interaction
Meadow_Trait_Variability_TM <- lmer(TM.MaxN~trait.co.var.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability_TM, test="Chisq") #no improvment without loss of seagrass variable
summary(Meadow_Trait_Variability_TM, correlation=FALSE)

Meadow_Trait_Variability_SST <- lmer(SST.MaxN~trait.co.var.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability_SST, test="Chisq") #Drop interaction
Meadow_Trait_Variability_SST <- lmer(SST.MaxN~trait.co.var.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability_SST, test="Chisq") #drop land use
Meadow_Trait_Variability_SST <- lmer(SST.MaxN~trait.co.var.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability_SST, test="Chisq")
summary(Meadow_Trait_Variability_SST, correlation=FALSE)

Meadow_Trait_Variability_LH <- lmer(LH.MaxN~trait.co.var.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability_LH, test="Chisq") #Drop interaction
Meadow_Trait_Variability_LH <- lmer(LH.MaxN~trait.co.var.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability_LH, test="Chisq") #no improvment without loss of seagrass variable
summary(Meadow_Trait_Variability_LH, correlation=FALSE)

Meadow_Trait_Variability_LV <- lmer(LV.MaxN~trait.co.var.pc1*land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability_LV, test="Chisq") #drop interaction
Meadow_Trait_Variability_LV <- lmer(LV.MaxN~trait.co.var.pc1+land.use.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability_LV, test="Chisq") #drop land use
Meadow_Trait_Variability_LV <- lmer(LV.MaxN~trait.co.var.pc1+Depth+(1|Site), data = Zanz_fish2, REML = FALSE)
drop1(Meadow_Trait_Variability_LV, test="Chisq")
summary(Meadow_Trait_Variability_LV, correlation=FALSE)

visreg(Meadow_Trait_Variability_LV)

#compare models
require(AICcmodavg)
TM_model_list<-list(Seagrass_Richness_TM, Seagrass_Cover_TM, Functional_Dispersion_TM, Functional_Richness_TM, Meadow_Traits_TM, Seagrass_Composition_TM, Meadow_Trait_Variability_TM)
model_names<-c("Seagrass Richness", "Seagrass Cover", "Functional Dispersion", "Functional Richness", "Meadow Traits", "Seagrass Composition", "Meadow Trait Variability")
modelsel<-aictab(TM_model_list, model_names, second.ord=T)
write.table(modelsel, "~/Dropbox/Phd/Paper II - Seagrass diversity fisheries/Data/TM_Model_Selection.txt", sep="\t")
modelsel
library(piecewiseSEM)
rsquared(TM_model_list)
anova(Seagrass_Richness_TM, Seagrass_Cover_TM, Functional_Dispersion_TM, Functional_Richness_TM, Meadow_Traits_TM, Seagrass_Composition_TM, Meadow_Trait_Variability_TM)
TM_model_list

#compare models
require(AICcmodavg)
SST_model_list<-list(Seagrass_Richness_SST, Seagrass_Cover_SST, Functional_Dispersion_SST, Functional_Richness_SST, Meadow_Traits_SST, Seagrass_Composition_SST, Meadow_Trait_Variability_SST)
model_names<-c("Seagrass Richness", "Seagrass Cover", "Functional Dispersion", "Functional Richness", "Meadow Traits", "Seagrass Composition", "Meadow Trait Variability")
modelsel<-aictab(SST_model_list, model_names, second.ord=T)
write.table(modelsel, "~/Dropbox/Phd/Paper II - Seagrass diversity fisheries/Data/SST_Model_Selection.txt", sep="\t")
modelsel
rsquared(SST_model_list)
SST_model_list

anova(Seagrass_Richness_SST, Seagrass_Cover_SST, Functional_Dispersion_SST, Functional_Richness_SST, Meadow_Traits_SST, Seagrass_Composition_SST, Meadow_Trait_Variability_SST)

#compare models
require(AICcmodavg)
LH_model_list<-list(Seagrass_Richness_LH, Seagrass_Cover_LH, Functional_Dispersion_LH, Functional_Richness_LH, Meadow_Traits_LH, Seagrass_Composition_LH, Meadow_Trait_Variability_LH)
model_names<-c("Seagrass Richness", "Seagrass Cover", "Functional Dispersion", "Functional Richness", "Meadow Traits", "Seagrass Composition", "Meadow Trait Variability")
modelsel<-aictab(LH_model_list, model_names, second.ord=T)
write.table(modelsel, "~/Dropbox/Phd/Paper II - Seagrass diversity fisheries/Data/LH_Model_Selection.txt", sep="\t")
modelsel
rsquared(LH_model_list)
LH_model_list

anova(Seagrass_Richness_LH, Seagrass_Cover_LH, Functional_Dispersion_LH, Functional_Richness_LH, Meadow_Traits_LH, Seagrass_Composition_LH, Meadow_Trait_Variability_LH)

#compare models
require(AICcmodavg)
LV_model_list<-list(Seagrass_Richness_LV, Seagrass_Cover_LV, Functional_Dispersion_LV, Functional_Richness_LV, Meadow_Traits_LV, Seagrass_Composition_LV, Meadow_Trait_Variability_LV)
model_names<-c("Seagrass Richness", "Seagrass Cover", "Functional Dispersion", "Functional Richness", "Meadow Traits", "Seagrass Composition", "Meadow Trait Variability")
modelsel<-aictab(LV_model_list, model_names, second.ord=T)
write.table(modelsel, "~/Dropbox/Phd/Paper II - Seagrass diversity fisheries/Data/LV_Model_Selection.txt", sep="\t")
modelsel
rsquared(LV_model_list)
LV_model_list

anova(Seagrass_Richness_LV, Seagrass_Cover_LV, Functional_Dispersion_LV, Functional_Richness_LV, Meadow_Traits_LV, Seagrass_Composition_LV, Meadow_Trait_Variability_LV)
#PLOTS

library(jtools)
library(broom)
library(ggstance)
library(broom.mixed)
plot_summs(Meadow_Trait_Variability_LV, Functional_Dispersion_LH, Meadow_Traits_SST, Meadow_Traits_TM, 
           scale = TRUE, inner_ci_level = .9, plot.distributions = TRUE,
           model.names = c("Low Value", "Local Household", "Small-Scale Trader", "Town Market"), rescale.distributions = FALSE, robust = TRUE, colors = "Dark2") + 
  theme(text = element_text(size = 14), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
  

display.brewer.all(colorblindFriendly = TRUE)

attach(Zanz_fish2)
ggplot(Zanz_fish2, aes(x=land.use.pc1, y=MaxN))+
  geom_point(aes(color=Site))+
  geom_smooth(method=lm)+
  scale_color_viridis(discrete = TRUE)

ggplot(Zanz_fish2, aes(x=land.use.pc1, y=Fishing.pressure)) + 
  geom_point()+
  geom_smooth(method=lm)+
  scale_color_viridis()

#nMDS
library(vegan)
community_matrix <- read.csv("Site x Species Matrix.csv", header = TRUE, sep=",")
community_trophic <- read.csv("Site x Species Matrix Trophic.csv", header = TRUE, sep=",")

community_matrix <- community_matrix[c(2:65)]

community_matrix<- community_matrix %>%
  group_by(site) %>% 
  summarise_all(mean)


library(textshape)
community_matrix <- column_to_rownames(community_matrix, 'site')

community.environment <-Zanz_fish2

community.environment<- community.environment %>%
  group_by(Site) %>% 
  summarise_all(mean)

community.environment <- column_to_rownames(community.environment, 'Site')
community.environment<- data.frame(community.environment[c(6,7,16,18)])


zanzibar_nmds <- metaMDS(community_matrix, distance = "bray")
zanzibar.envfit <- envfit(zanzibar_nmds, community.environment, permutations = 999)


site.scrs <- as.data.frame(scores(zanzibar_nmds, display = "species"))
site.scrs <- cbind(site.scrs, Family = rownames(site.scrs)) #add family names as variable if you want to display on plot
site.scrs <- cbind(site.scrs, Trophic.Group = community_trophic$Trophic.Feeder)
head(site.scrs)

env.scores.zanzibar <- as.data.frame(scores(zanzibar.envfit, display = "vectors")) #extracts relevant scores from envifit
env.scores.zanzibar <- cbind(env.scores.zanzibar, env.variables = rownames(env.scores.zanzibar)) #and then gives them their names

env.scores.zanzibar <- cbind(env.scores.zanzibar, pval = zanzibar.envfit$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores.zanzibar, pval<=0.05) #subset data to show variables significant at 0.05

head(env.scores.zanzibar)

#PLOT
library(viridis)
nmds.plot.zanzibar <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(site.scrs$Family), shape = factor(site.scrs$Trophic.Group)), size = 6) + #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Family", shape = "Trophic Group")+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10))+theme_grey()+ # add legend at right of plot
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) 


nmds.plot.zanzibar + labs(title = "Basic ordination plot") #displays plot

nmds.plot.zanzibar+
  geom_segment(data = env.scores.zanzibar, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant env variables
  ggrepel::geom_text_repel(data = env.scores.zanzibar, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 4, direction = "both", segment.size = 0.25)+ #add labels for env variables
  labs(title="Ordination with environmental vectors")

# function for ellipsess 
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the management factor
df_ell<- data.frame() #sets up a data frame before running the function.
for(g in levels(site.scrs$Trophic.Group)){
  df_ell<- rbind(df_ell, cbind(as.data.frame(with(site.scrs[site.scrs$Trophic.Group==g,],
                                                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,Trophic.Group=g))
}

df_ell <- data.frame()

for(g in levels(site.scrs$Trophic.Group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(site.scrs[site.scrs$Trophic.Group==g,],
                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),
                                                                                                    length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                ,Trophic.Group=g))
}

# data for labelling the ellipse
NMDS.mean.dune=aggregate(site.scrs[ ,c("NMDS1", "NMDS2")], 
                         list(group = site.scrs$Trophic.Group), mean)

# data for labelling the ellipse
NMDS.mean=aggregate(site.scrs[,c("NMDS1", "NMDS2")], 
                    list(group = site.scrs$Trophic.Group), mean)



nmds.plot.zanzibar+ 
  geom_path(data = df_ell.dune.family, aes(x = NMDS1, y = NMDS2)) #this is the ellipse, seperate ones by Site. 










data.scores = as.data.frame(scores(zanzibar_nmds))



library(viridis)
install.packages("viridis")







require(plyr)
community_matrix<-cbind(community_matrix[c()], predict(zanzibar_nmds, newdata = community_matrix))

assumeblage<-ddply(community_matrix, .(sample), summarize,  comp.pc1=mean(NMDS1), comp.pc2=mean(NMDS1))


  
  