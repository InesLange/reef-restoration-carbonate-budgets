####
## data analysis for:
## Coral restoration can drive rapid recovery of reef carbonate budgets
##
## Lange ID, Razak TB, Perry CT, Maulana PG, Prasetya ME, Irwan, Lamont TAC 
## 2023. under review

#
# Author Ines Lange
#
# This code compares coral cover, rugosity and carbonate budgets as well as 
# genera contributions and coral colony abundances and sizes for degraded, 
# restored and healthy reef sites in the Spermonde Archipelago, Indonesia.
# Reefs were restored with the Mars Assisted Reef Restoration System (MARRS).
# Data were collected in September 2022 and May 2023 
# by Ines Lange (ReefBudget surveys, coral growth, skeletal density) 
# and Mars Sustainable Solutions (MSS) survey team (fish and urchin surveys)

## load packages

# organise and plot data
library(tidyverse)
library(plyr)
library(viridis)
library(patchwork) #multipanel figures

# mixed models
library(lme4)
library(lmerTest)
library(emmeans) #pairwise comparison

# PCA
library(vegan)
library(FactoMineR)
library(factoextra)
library(pairwiseAdonis)

#####
## load data
#####

# coral cover, rugosity, carbonate budget data
marrs <- read.csv("data/MARRS_ReefBudget.csv")%>%
  mutate(status = as.factor(status))
marrs$status <- factor(marrs$status, 
                       levels = c("degraded", "0", "1","2","4","healthy"))

# coral colony abundance and size data
colony <- read.csv("data/MARRS_colonysize.csv", strip.white = T) %>%
  mutate(status = as.factor(status))
colony$status <- factor(colony$status, levels = c("degraded", "0", "1","2","4", "healthy"))


#####
## summarise carbonate budget data at site level and by status
#####

avg_site <- marrs %>%
  group_by(status,site) %>%
  summarise_at(vars(9:20),list(mean=mean,sd=sd))
avg_site <- as.data.frame(avg_site)
write.csv(avg_site, "avg_site.csv", row.names=F)

avg_status <- avg_site %>%
  group_by(status) %>%
  summarise_at(vars(2:13),list(avg=mean,sd=sd))
avg_status <- as.data.frame(avg_status)
write.csv(avg_status, "avg_status.csv", row.names=F)

#####
## summarise colony abundance and size data at transect and site level
#####

avg_colony <- ddply(colony, c("status", "site", "transect"), summarise,
               abund    = as.numeric(length(size)),
               size = mean(size),
               sd   = sd(size),
               se   = sd / sqrt(abund),
               median = median(size),
               skew = 3*(size-median)/sd)
avg_colony <- as.data.frame(avg_colony)

avg_colony_site <- avg_colony %>%
  group_by(status,site) %>%
  summarise_at(vars(2:3),list(mean=mean,sd=sd))
avg_colony_site <- as.data.frame(avg_colony_site)


#####
## linear mixed models to compare metrics between reefs with different restoration status
#####

cover.lmer <- lmer(coral_cover ~  status + (1|site), data=marrs) #site as random
summary(cover.lmer) 
anova(cover.lmer) #status sign impacts coral cover F5,12=78.169, p<0.001
plot(cover.lmer) #ok
qqnorm(resid(cover.lmer))
qqline(resid(cover.lmer)) #ok
emmeans(cover.lmer, "status") # overlap CL 1-2, 4-healthy
pairs(emmeans(cover.lmer, "status")) # all diff except 0=1, 1=2, 4=hea

rugosity.lmer <- lmer(rugosity ~  status + (1|site), data=marrs) #site as random
summary(rugosity.lmer) 
anova(rugosity.lmer) #status sign impacts rugosity F5,12=15.966, p<0.001
plot(rugosity.lmer) #ok
qqnorm(resid(rugosity.lmer))
qqline(resid(rugosity.lmer)) #ok
emmeans(rugosity.lmer, "status") 
pairs(emmeans(rugosity.lmer, "status")) # all diff except deg=0, 0=1, 1=2=4=hea

gross_prod.lmer <- lmer(gross_prod ~  status + (1|site), data=marrs) #site as random
summary(gross_prod.lmer) 
anova(gross_prod.lmer) #status sign impacts gross_prod F5,12=15.966, p<0.001
plot(gross_prod.lmer) #ok
qqnorm(resid(gross_prod.lmer))
qqline(resid(gross_prod.lmer)) #ok
emmeans(gross_prod.lmer, "status") # overlap CL 0-1, 4-hea
pairs(emmeans(gross_prod.lmer, "status")) # all diff except 0=1, 4=hea

gross_ero.lmer <- lmer(gross_ero ~  status + (1|site), data=marrs) #site as random
summary(gross_ero.lmer) 
anova(gross_ero.lmer) #status does NOT sign impact gross_ero F5,12=1.683, p=0.213
plot(gross_ero.lmer) #ok
qqnorm(resid(gross_ero.lmer))
qqline(resid(gross_ero.lmer)) #ok

net_prod.lmer <- lmer(net_prod ~  status + (1|site), data=marrs) #site as random
summary(net_prod.lmer) 
anova(net_prod.lmer) #status sign impacts net_prod F5,12=69.478, p<0.001
plot(net_prod.lmer) #ok
qqnorm(resid(net_prod.lmer))
qqline(resid(net_prod.lmer)) #ok
emmeans(net_prod.lmer, "status") # overlap CL 0-1, 4-hea
pairs(emmeans(net_prod.lmer, "status")) # all diff except 0=1, 4=hea

RAPmax.lmer <- lmer(RAPmax ~  status + (1|site), data=marrs) #site as random
summary(RAPmax.lmer) 
anova(RAPmax.lmer) #status sign impacts RAPmax F5,12=44.734, p<0.001
plot(RAPmax.lmer) #ok
qqnorm(resid(RAPmax.lmer))
qqline(resid(RAPmax.lmer)) #ok
emmeans(RAPmax.lmer, "status") # overlap CL 0-1, 2-4, 1-hea, 2-hea, 4=hea
pairs(emmeans(RAPmax.lmer, "status")) # all diff except 0=1, 1=hea, 2=hea, 2=4, 4=hea 

#####
## Figure 1
## plot recovery of metrics over time (transect-level data: boxplots, site means: jitter points) 
#####

cover.plot <- ggplot(marrs, aes(x = status, y = coral_cover, fill=status))
(p1 <- cover.plot + geom_boxplot(show.legend=FALSE) + 
  geom_jitter(data=avg_site, aes(x = status, y = coral_cover_mean), 
              size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
  ylab("Coral cover (%)") + xlab("") + theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12)))

rugosity.plot <- ggplot(marrs, aes(x = status, y = rugosity, fill=status))
(p2 <- rugosity.plot + geom_boxplot(show.legend=FALSE) + 
    geom_jitter(data=avg_site, aes(x = status, y = rugosity_mean), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
    ylab("Rugosity") + xlab("") + theme_classic() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)))

abund.plot <- ggplot(avg_colony, aes(x = status, y = abund, fill=status))
(p3 <- abund.plot + geom_boxplot(show.legend=FALSE) + 
    geom_jitter(data=avg_colony_site, aes(x = status, y = abund_mean), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
    ylab("Colony abund. (ind. transect-1)") + xlab("") + theme_classic() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)))

size.plot <- ggplot(avg_colony, aes(x = status, y = size, fill=status))
(p4 <- size.plot + geom_boxplot(show.legend=FALSE) + 
    geom_jitter(data=avg_colony_site, aes(x = status, y = size_mean), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
    ylab("Colony size (cm)") + xlab("") + theme_classic() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)))

gross_prod.plot <- ggplot(marrs, aes(x = status, y = gross_prod, fill=status))
(p5 <- gross_prod.plot + geom_boxplot(show.legend=FALSE) + 
    geom_jitter(data=avg_site, aes(x = status, y = gross_prod_mean), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
    ylab("Gross production (kg m-2 yr-1)") + xlab("") + theme_classic() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)))

gross_ero.plot <- ggplot(marrs, aes(x = status, y = gross_ero, fill=status))
(p6 <- gross_ero.plot + geom_boxplot(show.legend=FALSE) + 
    geom_jitter(data=avg_site, aes(x = status, y = gross_ero_mean), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
    ylab("Gross erosion (kg m-2 yr-1)") + xlab("") + theme_classic() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)))

net_prod.plot <- ggplot(marrs, aes(x = status, y = net_prod, fill=status))
(p7 <- net_prod.plot + geom_boxplot(show.legend=FALSE) + 
    geom_jitter(data=avg_site, aes(x = status, y = net_prod_mean), 
                size=3, colour="black", fill="white", pch=21, width=0.1, show.legend=FALSE) +
    ylab("Net budget (kg m-2 yr-1)") + xlab("Time since restoration (years)") + theme_classic() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)))

# combine plots into Fig 1
Fig1_recovery <- p1+p2+p3+p4+p5+p6+p7 + plot_layout(ncol=2) +
  plot_annotation(tag_levels="A")
Fig1_recovery

ggsave("figures/Fig1_recovery.pdf", width = 10, height = 12)


#####
## Analyse contribution to coral carbonate production by different coral genera
#####

#transform genus-level contributions to long format
genus_long <- marrs[,c(1:3,7,23:37)] %>%
  pivot_longer(cols=c("ACRA", "ACRO", "ACRT", "POCB", "STYB", "SERB", "PORB", 
                      "HCB", "PORM", "ISOS", "HCM", "HCE", "HCP", "FUN", "Other"), 
               names_to="genus", values_to="prod")

#calculate averages for each status
genus_avg <- genus_long %>%
  group_by(status,genus) %>%
  summarise_at(vars(4),list(prod=mean,sd=sd))
genus_avg <- as.data.frame(genus_avg)

genus_order <- c("ACRA", "ACRO", "ACRT",  "HCB", "POCB", "STYB", "SERB", "PORB", 
                 "PORM", "ISOS", "HCM", "HCE", "HCP", "FUN", "Other")
genus_avg$genus <- factor(genus_avg$genus, levels = rev(genus_order))

# plot without degraded sites as very little coral carbonate production
(p8 <- ggplot(genus_avg[(!(genus_avg$status=="degraded")), ], aes(x = status, y = prod, fill = genus)) +
    scale_fill_viridis(option="C", discrete=TRUE)+
    geom_bar(stat = "identity", position="fill") + #position="fill" to display 100% for each
    labs(x = "Time since restoration (years)", y = "Carbonate production (%)") +
    theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=12), 
                       legend.title = element_text(size=12),legend.text = element_text(size=12)))


#####
## PCA for genera-level contributions to total coral carbonate production
#####

coralPCA <- marrs[-7,] #remove one transect without any coral cover (as causing error in PERMANOVA)
coralPCA2 <- coralPCA[(!(coralPCA$status=="degraded")), ] # run PCA without degraded

marrs.spe_coral <- coralPCA2[,23:37] # specify taxa columns
marrs.spe_coral[marrs.spe_coral < 0.001] <- 0
marrs.spe_normcoral <- as.data.frame(sqrt(marrs.spe_coral[,])) # normalize data
#marrs.spe_normcoral[is.na(marrs.spe_normcoral)] <- 0 #replace NA with 0
marrs.env <- subset(coralPCA2[7]) # status as grouping variable

#Run PCA
marrs.pca_coral<-rda(marrs.spe_normcoral, scale=TRUE)
marrs.pca_coral
summary(marrs.pca_coral)

#how many PC
screeplot(marrs.pca_coral, type="lines") #PC1 and PC2 explain most variability

# which coral genera are driving differences among sites?
ef.pca<-envfit(marrs.pca_coral,marrs.spe_coral,permu=999) 
ef.pca #all genera except POCB and PORB sign drive differences among sites (p<0.05*)

# PERMANOVA to assess differences between sites with different restoration status
permanova <- adonis2(marrs.spe_coral~status, data=coralPCA2)
permanova # time since restoration significantly affects genera contributions to Coral G
posthoc_coral <- pairwise.adonis(marrs.spe_coral, factors=coralPCA2$status)
posthoc_coral # all diff except 2=4 (0=1 and 4=hea)

#Permdisp
permdisp <- betadisper(vegdist(marrs.spe_normcoral, method="bray", na.rm=TRUE), coralPCA2$status) 
par(mfrow=c(1,2))
plot(permdisp, main="PCoA")
boxplot(permdisp, main="Distance to centroids")
anova(permdisp) # ns, meaning group dispersion is homogenous 
# this means differences between centroids are solely due to differences in composition

# plot PCA with status as grouping
marrs.pca_coral <- PCA(marrs.spe_normcoral, graph=FALSE)
p9 <- fviz_pca_biplot(marrs.pca_coral,
                          habillage = coralPCA2$status,
                          addEllipses = TRUE, ellipse.level = 0.9,
                          label = "var", col.var = "black", repel = TRUE,
                          #select.var = list(contrib=10) # only display 10 most contributing taxa
                          title="")
(p9 <- p9 + scale_color_manual(values=c("#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) + 
    theme_classic() + theme(legend.position="top") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))) 

# combine stacked bars and PCA into Figure 2 
Fig2_coral_prod <- p8+p9+
  plot_annotation(tag_levels="A")
Fig2_coral_prod

ggsave("figures/Fig2_coral_prod.pdf", width = 10, height = 4.7)

