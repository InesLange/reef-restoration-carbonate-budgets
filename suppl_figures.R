####
## Supplementary Figures for:
## Coral restoration can drive rapid recovery of reef carbonate budgets
##
## Lange ID, Razak TB, Perry CT, Maulana PG, Prasetya ME, Irwan, Lamont TAC 
## Current Biology (2024)

#
# Author Ines Lange

# load packages
library(tidyverse)
library(patchwork) #multipanel figures

# PCA
library(vegan)
library(FactoMineR)
library(factoextra)
library(pairwiseAdonis)

# mixed models
library(lme4)
library(lmerTest)
library(emmeans) #pairwise comparison

#evenness
library(adiv)

#####
## Figure S1A: PCA benthic community composition
#####

# load data
benthicPCA <- read.csv("data/MARRS_benthic.csv")%>%
  mutate(status = as.factor(status))
benthicPCA$status <- factor(benthicPCA$status, levels = c("degraded", "0", "1","2","4", "healthy"))

marrs.spe_benthic <- benthicPCA[,5:26] # specify taxa columns
marrs.spe_benthic[marrs.spe_benthic < 0.001] <- 0
marrs.spe_normbenthic <- as.data.frame(sqrt(marrs.spe_benthic[,])) # normalize data
#marrs.spe_normcoral[is.na(marrs.spe_normcoral)] <- 0 #replace NA with 0
marrs.env <- subset(benthicPCA[4]) # status as grouping variable

#Run PCA
marrs.pca_benthic<-rda(marrs.spe_normbenthic, scale=TRUE)
marrs.pca_benthic
summary(marrs.pca_benthic)

#how many PC
screeplot(marrs.pca_benthic, type="lines") #PC1 explains most variability

# which benthic genera are driving differences among sites?
ef.pca<-envfit(marrs.pca_benthic,marrs.spe_benthic,permu=999) 
ef.pca #all benthic groups except SOC, MA and HCE sign drive differences among sites (p<0.05)

# PERMANOVA to assess differences between sites with different restoration status
permanova2 <- adonis2(marrs.spe_benthic~status, data=benthicPCA)
permanova2 # time since restoration significantly affects benthic community composition
posthoc_benthic <- pairwise.adonis(marrs.spe_benthic, factors=benthicPCA$status)
posthoc_benthic # all different except 0=1 (p<0.05 but >0.01)

#Permdisp
permdisp <- betadisper(vegdist(marrs.spe_normbenthic, method="bray", na.rm=TRUE), benthicPCA$status) 
par(mfrow=c(1,2))
plot(permdisp, main="PCoA")
boxplot(permdisp, main="Distance to centroids")
anova(permdisp) # sign, dispersion is higher among mature restoration and healthy reef sites

# plot PCA with status as grouping
marrs.pca_benthic <- PCA(marrs.spe_normbenthic, graph=FALSE)
p13 <- fviz_pca_biplot(marrs.pca_benthic,
                       habillage = benthicPCA$status,
                       addEllipses = TRUE, ellipse.level = 0.9,
                       label = "var", col.var = "black", repel = TRUE,
                       #select.var = list(contrib=10) # only display 10 most contributing taxa
                       title="")
(p13 <- p13 +  
    theme_classic() + theme(legend.position="right") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))) 

#####
## Fig S1 B&C: Shannon-Wiener-Diversity Index and Smith-Wilson Evenness index of coral genera composition
#####

#calculate Shannon-Wiener-Diversity Index and Smith-Wilson Evenness index
cover_diversity <- benthicPCA[,12:26]
cover_diversity[cover_diversity < 0] <- 0 

cover_shannon <- data.frame(diversity(cover_diversity, index="shannon"))
cover_shannon <- cbind(benthicPCA[,1:4],cover_shannon)
cover_shannon$shannon <- cover_shannon[,5] #change name of column to shannon
cover_shannon <- cover_shannon[!(cover_shannon$status %in% "degraded"),] #remove degraded sites as hardly any coral

cover_smithwilson <- data.frame(specieseve(cover_diversity)) # 2 transects without coral removed from analysis
cover_smithwilson <- cbind(benthicPCA[-c(7,8),1:4],cover_smithwilson) #remove 2 transects without coral from benthic PCA
cover_smithwilson <- cover_smithwilson[!(cover_smithwilson$status %in% "degraded"),] #remove degraded sites as hardly any coral

# plot differences across sites
cover_shannon_mean <- cover_shannon %>% 
  group_by(island, site, status) %>% 
  summarise_at(vars(2),list(avg=mean,sd=sd))
cover_shannon_mean <- as.data.frame(cover_shannon_mean)

p14 <- ggplot(cover_shannon, aes(x = status, y = shannon, fill=status))
(p14 <- p14 + geom_boxplot(show.legend=FALSE) + 
        scale_fill_manual(values=c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
    ylab("Shannon-Wiener Diversity Index") + xlab("Time since restoration (years)") + theme_classic() + 
    geom_jitter(data=cover_shannon_mean, aes(x = status, y = avg, fill=status), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
        theme(axis.text=element_text(size=12), axis.title=element_text(size=14)))
  
cover_smithwilson_mean <- cover_smithwilson %>% 
  group_by(island, site, status) %>% 
  summarise_at(vars(7),list(avg=mean,sd=sd))
cover_smithwilson_mean <- as.data.frame(cover_smithwilson_mean)

p15 <- ggplot(cover_smithwilson, aes(x = status, y = SmithWilson, fill=status))
(p15 <- p15 + geom_boxplot(show.legend=FALSE) + 
    scale_fill_manual(values=c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
    ylab("Smith-Wilson Evenness Index") + xlab("Time since restoration (years)") + theme_classic() + 
    geom_jitter(data=cover_smithwilson_mean, aes(x = status, y = avg, fill=status), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14)))

# test differences between sites
shannon.lmer <- lmer(shannon ~  status + (1|site), data=cover_shannon) #site as random
summary(shannon.lmer) #isSingular - overfitted?, ANOVA at site-level instead
summary(aov(avg~status, data=cover_shannon_mean)) # F4,10=7.593, p=0.004
TukeyHSD(aov(avg~status, data=cover_shannon_mean))
multcompLetters4(aov(avg~status, data=cover_shannon_mean), TukeyHSD(aov(avg~status, data=cover_shannon_mean)), reversed=TRUE)
#a,b,ab,b,b

smithwilson.lmer <- lmer(SmithWilson ~  status + (1|site), data=cover_smithwilson) #site as random
summary(smithwilson.lmer) #isSingular - overfitted?, ANOVA at site-level instead
summary(aov(avg~status, data=cover_smithwilson_mean)) # F4,10=5.521, p=0.013
TukeyHSD(aov(avg~status, data=cover_smithwilson_mean))
multcompLetters4(aov(avg~status, data=cover_smithwilson_mean), TukeyHSD(aov(avg~status, data=cover_smithwilson_mean)), reversed=FALSE)
# ab,a,b,b,b

#combine plots into Fig S2
FigS2_diversity <- (p13+plot_spacer())/(p14 + p15) +
  plot_annotation(tag_levels="a")
FigS2_diversity

ggsave("figures/FigS2_diversity.pdf", width = 10, height = 4)


#####
## Fig S3: Mean and median colony size and skewness and kurtosis of size distributions
#####

# plot size distribution data of all coral colonies

# mean colony size
size.plot <- ggplot(avg_colony[!(avg_colony$mean_size=="0"),], aes(x = status, y = mean_size, fill=status))
(p4 <- size.plot + geom_boxplot(show.legend=FALSE) + 
    geom_jitter(data=avg_colony_site, aes(x = status, y = mean_size_mean), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
    ylab("Mean colony size (cm)") + xlab("") + theme_classic() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)))
  
# median colony size
median.plot <- ggplot(avg_colony, aes(x = status, y = median, fill=status))
(p10 <- median.plot + geom_boxplot(show.legend=FALSE) + 
    geom_jitter(data=avg_colony_site, aes(x = status, y = median_mean), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
    ylab("Median colony size (cm)") + xlab("") + theme_classic() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)))

# standard deviation colony size
SD.plot <- ggplot(avg_colony[!(avg_colony$mean_size=="0"),], aes(x = status, y = SD, fill=status))
(p0 <- SD.plot + geom_boxplot(show.legend=FALSE) + 
    geom_jitter(data=avg_colony_site, aes(x = status, y = SD_mean), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
    ylab("SD colony size (cm)") + xlab("") + theme_classic() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)))

# skewness
skew.plot <- ggplot(avg_colony[!(avg_colony$mean_size=="0"),], aes(x = status, y = skew, fill=status))
(p11 <- skew.plot + geom_boxplot(show.legend=FALSE) + 
    geom_jitter(data=avg_colony_site, aes(x = status, y = skew_mean), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
    ylab("Skew of size distribution") + xlab("") + theme_classic() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)))

# kurtosis
kurtosis.plot <- ggplot(avg_colony[!(avg_colony$mean_size=="0"),], aes(x = status, y = kurtosis, fill=status))
(p12 <- kurtosis.plot + geom_boxplot(show.legend=FALSE) + 
    geom_jitter(data=avg_colony_site, aes(x = status, y = kurtosis_mean), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
    ylab("kurtosis of size distribution") + xlab("") + theme_classic() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)))

# combine plots into Figure S3 
FigS2_colonySize <- p4+p10+p11+p12 + plot_layout(ncol=2) +
  plot_annotation(tag_levels="a")
FigS2_colonySize

ggsave("figures/FigS2_colonySize.pdf", width = 10, height = 7)

###
# test differences
###

mean_size.lmer <- lmer(mean_size ~  status + (1|site), data=avg_colony[!(avg_colony$mean_size=="0"),]) #remove transect with 0 colonies
summary(mean_size.lmer) 
anova(mean_size.lmer) #status sign impacts coral mean_size F5,12=24.613, p<0.001
plot(mean_size.lmer) #ok
qqnorm(resid(mean_size.lmer))
qqline(resid(mean_size.lmer)) #ok
emmeans(mean_size.lmer, "status") # overlap CL degr-0-1, 1-2, 2-4, 4-healthy
pairs(emmeans(mean_size.lmer, "status")) # all diff except degr=0=1, degr=1=2, 2=4, 4=hea

median_size.lmer <- lmer(median ~  status + (1|site), data=avg_colony[!(avg_colony$median=="0"),]) #remove transect with 0 colonies
summary(median_size.lmer) 
anova(median_size.lmer) #status sign impacts coral median_size F5,12=12.111, p<0.001
plot(median_size.lmer) #ok
qqnorm(resid(median_size.lmer))
qqline(resid(median_size.lmer)) #ok
emmeans(median_size.lmer, "status") # overlap CL degr-0-1, 1-2, 2-4-healthy
pairs(emmeans(median_size.lmer, "status")) # all diff except degr=0=1, 1=2=4=hea

SD.lmer <- lmer(SD ~  status + (1|site), data=avg_colony[!(avg_colony$mean_size=="0"),]) #site as random
summary(SD.lmer) 
anova(SD.lmer) #status sign impacts SD of mean colony size F5,12=18.795, p<0.001
plot(SD.lmer) #ok
qqnorm(resid(SD.lmer))
qqline(resid(SD.lmer)) #ok
emmeans(SD.lmer, "status") # overlap CL deg-0-1, 1-2, 2-4, 4-healthy
pairs(emmeans(SD.lmer, "status")) # all diff except degr=0=1=2,2=4, 4=hea

skew.lmer <- lmer(skew ~ status + (1|site), data=avg_colony[complete.cases(avg_colony), ]) #remove transects with error
summary(skew.lmer) #isSingular - overfitted?, ANOVA at site-level instead
summary(aov(skew_mean~status, data=avg_colony_site)) # F5,11=8.809, p=0.001
TukeyHSD(aov(skew_mean~status, data=avg_colony_site))
multcompLetters4(aov(skew_mean~status, data=avg_colony_site), TukeyHSD(aov(skew_mean~status, data=avg_colony_site)), reversed=TRUE)
# a,ab,bc,ab,abc,c

kurtosis.lmer <- lmer(kurtosis ~ status + (1|site), data=avg_colony[complete.cases(avg_colony), ]) #remove transects with error
summary(kurtosis.lmer) #isSingular - overfitted? ANOVA at site-level instead
summary(aov(kurtosis_mean~status, data=avg_colony_site)) # F5,11=9.402, p=0.001
TukeyHSD(aov(kurtosis_mean~status, data=avg_colony_site))
multcompLetters4(aov(kurtosis_mean~status, data=avg_colony_site), TukeyHSD(aov(kurtosis_mean~status, data=avg_colony_site)), reversed=TRUE)
# a,ab,bc,ab,abc,c


#####
## Figure S3: Coral linear growth and skeletal density
#####

#load data
growth <- read.csv("data/MARRS_coralgrowth.csv", strip.white = T) %>%
  mutate(taxa = as.factor(taxa))
growth <- growth[(!(growth$taxa=="ACRH"| growth$taxa=="SER")),] # remove as only 1 replicate/taxa
growth$taxa <- factor(growth$taxa, levels = c("ACRA", "ACRO/T", "ACRO","STY","POC", "PORB", "HYD", "ISO"))

p.growth <- ggplot(growth, aes(x = taxa, y = growth))
(p.growth <- p.growth + geom_boxplot() + ylab("Linear growth (cm/yr)") + xlab("Genera") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) + theme_classic())

p.dens <- ggplot(growth, aes(x = taxa, y = density))
(p.dens <- dens + geom_boxplot() + ylab("Skeletal density (g/cm3)")+ xlab("Genera") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) + theme_classic())

# correlation between growth and density
growth_density <- ggplot(growth, aes(x = growth, y = density))
growth_density + geom_point() +  geom_smooth(method=lm, se=TRUE) + 
  ylab("Skeletal density (g/cm3)") + xlab("Linear growth (cm/yr)") + theme_classic()

#combine growth and density plots
FigS3 <- p.growth + p.dens +
  plot_annotation(tag_levels="a")
FigS3

#save plot
ggsave("figures/FigS3_coralgrowth.pdf", width = 8, height = 3.5)

