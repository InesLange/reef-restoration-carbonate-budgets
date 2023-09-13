####
## Supplementary Figures for:
## Coral restoration can drive rapid recovery of reef carbonate budgets
##
## Lange ID, Razak TB, Perry CT, Maulana PG, Prasetya ME, Irwan, Lamont TAC 
## 2023. under review

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

#####
## Figure S2: PCA benthic community composition
#####

# load data
coverPCA <- read.csv("data/MARRS_benthic.csv")%>%
  mutate(status = as.factor(status))
coverPCA$status <- factor(coverPCA$status, levels = c("degraded", "0", "1","2","4", "healthy"))

marrs.spe_cover <- coverPCA[,5:26] # specify taxa columns
marrs.spe_cover[marrs.spe_cover < 0.001] <- 0
marrs.spe_normcover <- as.data.frame(sqrt(marrs.spe_cover[,])) # normalize data
#marrs.spe_normcover[is.na(marrs.spe_normcover)] <- 0 #replace NA with 0
marrs.env_cover <- subset(coverPCA[4]) # status as grouping variable

#Run PCA
marrs.pca_cover<-rda(marrs.spe_normcover, scale=TRUE)
marrs.pca_cover
summary(marrs.pca_cover)

#how many PC
screeplot(marrs.pca_cover, type="lines") #PC1 and PC2 explain most variability

# which cover genera are driving differences among sites?
ef.pca2<-envfit(marrs.pca_cover,marrs.spe_cover,permu=999) 
ef.pca2 #all benthic groups except SOC, MA and HCE sign drive differences among sites (p<0.05*)

# PERMANOVA to assess differences between sites with different restoration status
permanova2 <- adonis2(marrs.spe_cover~status, data=coverPCA)
permanova2 # time since restoration significantly affects benthi community composition
posthoc_cover <- pairwise.adonis(marrs.spe_cover, factors=coverPCA$status)
posthoc_cover # all diff except 0=1

#Permdisp
permdisp2 <- betadisper(vegdist(marrs.spe_normcover, method="bray", na.rm=TRUE), coverPCA$status) 
par(mfrow=c(1,2))
plot(permdisp2, main="PCoA")
boxplot(permdisp2, main="Distance to centroids")
anova(permdisp2) # sign, meaning group dispersion is not homogenous 
TukeyHSD(permdisp2) # 4> deg,1,2, hea> deg 
# this means beta diversity may be higher at mature restored and healthy sites

# plot PCA with status as grouping
marrs.pca_cover <- PCA(marrs.spe_normcover, graph=FALSE)
p10 <- fviz_pca_biplot(marrs.pca_cover,
                      habillage = coverPCA$status,
                      addEllipses = TRUE, ellipse.level = 0.9,
                      label = "var", col.var = "black", repel = TRUE,
                      #select.var = list(contrib=10) # only display 10 most contributing taxa
                      title="")
(p10 <- p10 + scale_color_manual(values=c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) + 
    theme_classic() + theme(legend.position="top") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))) 

#save plot
ggsave("figures/FigS2_PCAcover.pdf", width = 5, height = 4.5)


#####
## Figure S3: Size-frequency distributions of coral colonies
#####

# load data
colony <- read.csv("data/MARRS_colonysize.csv", strip.white = T) %>%
  mutate(status = as.factor(status))
colony$status <- factor(colony$status, levels = c("degraded", "0", "1","2","4", "healthy"))

# adjust and group morphotypes
colony2 <- colony %>%
  mutate(morph2 = gsub("corymbose", "corymbose Acropora", morphology))
colony2 <- colony2 %>%
  mutate(morph2 = gsub("digitate", "corymbose Acropora", morph2))
colony2 <- colony2 %>%
  mutate(morph2 = gsub("table", "corymbose Acropora", morph2))
colony2 <- colony2 %>%
  mutate(morph2 = gsub("arborescent", "arborescent Acropora", morph2))
colony2 <- colony2 %>%
  mutate(morph2 = gsub("caespitose", "arborescent Acropora", morph2))
colony2 <- colony2 %>%
  mutate(morph2 = gsub("hispidose", "arborescent Acropora", morph2))
colony2 <- colony2 %>%
  mutate(morph2 = gsub("columnar", "massive", morph2))
colony2 <- colony2 %>%
  mutate(morph2 = gsub("submassive", "massive", morph2))
colony2 <- colony2 %>%
  mutate(morph2 = gsub("freeliving", "other", morph2))
colony2 <- colony2 %>%
  mutate(morph2 = gsub("plating", "other", morph2))
colony2 <- colony2 %>%
  mutate(morph2 = gsub("foliose", "other", morph2))

colony2$morph2 <- as.factor(colony2$morph2)

# plot Size-Frequency distributions
# one line for each group as different number of sites

pal_freq <- c("#F8766D","#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")

(p11 <- ggplot(data=colony2, aes(size, y = after_stat(count)/12, colour=status, fill=status)) +
    geom_density(alpha = .5) +
    geom_line(stat = "density") +
    scale_x_continuous(trans = "log10", limits = c(0.9,300)) + 
    scale_y_continuous(limits = c(0,30), position="left") +
    facet_wrap(factor(morph2, levels=c('arborescent Acropora', 'corymbose Acropora', 'branching', 'massive', 'encrusting', 'other'))~., 
               ncol = 3, strip.position = "top") + #, scales = "free_y"
    theme_minimal() + theme(legend.position="top") +
    theme(text = element_text(colour = "black", size=14),
          panel.grid.minor = element_blank(),
          axis.line = element_line(),
          axis.ticks = element_line()) +
    labs(x = "log colony size (cm)", y = "Colonies per transect")) #
  

# save plot
ggsave("figures/FigS3_SizeFreq.pdf", width = 10, height = 7)


#####
## Figure S4: Coral linear growth and skeletal density
#####

#load data
growth <- read.csv("data/MARRS_coralgrowth.csv", strip.white = T) %>%
  mutate(taxa = as.factor(taxa))
growth <- growth[(!(growth$taxa=="ACRH"| growth$taxa=="SER")),] # remove as only 1 replicate/taxa
growth$taxa <- factor(growth$taxa, levels = c("ACRA", "ACRO/T", "ACRO","STY","POC", "PORB", "HYD", "ISO"))

p.dens <- ggplot(growth, aes(x = taxa, y = density))
(p.dens <- dens + geom_boxplot() + ylab("Skeletal density (g/cm3)")+ xlab("Genera") + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) + theme_classic())

p.growth <- ggplot(growth, aes(x = taxa, y = growth))
(p.growth <- p.growth + geom_boxplot() + ylab("Linear growth (cm/yr)") + xlab("Genera") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) + theme_classic())

# correlation between growth and density
growth_density <- ggplot(growth, aes(x = growth, y = density))
growth_density + geom_point() +  geom_smooth(method=lm, se=TRUE) + 
  ylab("Skeletal density (g/cm3)") + xlab("Linear growth (cm/yr)") + theme_classic()

#combine growth and density plots
FigS4 <- p.growth + p.dens +
  plot_annotation(tag_levels="a")
FigS4

#save plot
ggsave("figures/FigS4_coralgrowth.pdf", width = 8, height = 3.5)
