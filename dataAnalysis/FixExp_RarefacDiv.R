
############################################################################################
    ###  Fixative Experiment - amplicon analysis  ###
############################################################################################

# This script: calculate rarefaction and diversity indices 

library(iNEXT)
library(phyloseq)
library(olsrr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(gtools)


############################################################################################
   ###  RAREFACTION AND COVERAGE  ###
############################################################################################

iNEXT.bac <- otu_table(
  bac.otu, taxa_are_rows=F)
iNEXT.euk <- otu_table(
  euk.otu, taxa_are_rows=F)

iNEXT.bac <- iNEXT(
  as.data.frame(otu_table(iNEXT.bac)), q=c(0),
  datatype="abundance", conf = 0.95, nboot = 100)
iNEXT.euk <- iNEXT(
  as.data.frame(otu_table(iNEXT.euk)), q=c(0),
  datatype="abundance", conf = 0.95, nboot = 100)

###################################

## RAREFACTION ##

rarefac.bac <- fortify(iNEXT.bac, type=1) %>%
  left_join(ENV, by=c(
    "site"="internal_id"))

rarefac.point.bac <- rarefac.bac[which(
  rarefac.bac$method == "observed"),]
rarefac.line.bac <- rarefac.bac[which(
  rarefac.bac$method != "observed"),]
rarefac.line.bac$method <- factor(rarefac.line.bac$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

rarefac.euk <- fortify(iNEXT.euk, type=1) %>%
  left_join(ENV, by=c(
    "site"="internal_id"))

rarefac.point.euk <- rarefac.euk[which(
  rarefac.euk$method == "observed"),]
rarefac.line.euk <- rarefac.euk[which(
  rarefac.euk$method != "observed"),]
rarefac.line.euk$method <- factor(rarefac.line.euk$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

###################################

## COVERAGE ##

cover.bac <- fortify(iNEXT.bac, type=2)  %>%
  left_join(ENV, by=c(
    "site"="internal_id"))
cover.euk <- fortify(iNEXT.euk, type=2) %>%
  left_join(ENV, by=c(
    "site"="internal_id"))

cover.point.bac <- cover.bac [which(
  cover.bac$method == "observed"),]
cover.line.bac <- cover.bac [which(
  cover.bac$method != "observed"),]
cover.line.bac$method <- factor(cover.line.bac$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

cover.point.euk <- cover.euk [which(
  cover.euk$method == "observed"),]
cover.line.euk <- cover.euk [which(
  cover.euk$method != "observed"),]
cover.line.euk$method <- factor(cover.line.euk$method,
    c("interpolated","extrapolated"),
    c("interpolation","extrapolation"))

coverage.bac <- ggplot(cover.bac, 
  aes(x=x, y=y, colour=site))+ 
  geom_line(aes(linetype = method), 
  lwd = 0.5, data = cover.line.bac) +
  scale_colour_discrete(guide = F) +
  scale_x_continuous(
    limits = c(0,1e+5)) +
  scale_y_continuous(
    breaks = seq(0.9,1,0.05), 
    limits = c(0.9,1)) +
  facet_grid(~treatment)+ 
  labs(x="Sample size", y="Sample coverage") +
  theme_bw(base_size = 12) + 
  theme(legend.position="none",
        axis.ticks = element_blank())

coverage.euk <- ggplot(cover.euk, 
  aes(x=x, y=y, colour=site))+ 
  geom_line(aes(linetype = method), 
  lwd = 0.5, data = cover.line.euk) +
  scale_colour_discrete(guide = F) +
  scale_x_continuous(
    limits = c(0,1e+5)) +
  scale_y_continuous(
    breaks = seq(0.9,1,0.05), 
    limits = c(0.9,1)) +
  facet_grid(~treatment)+ 
  labs(x="Sample size", y="Sample coverage") +
  theme_bw(base_size = 12) + 
  theme(legend.position="none",
        axis.ticks = element_blank())

###################################

## RICHNESS

richness.bac <- ggplot(rarefac.bac, 
  aes(x=x, y=y, colour=site)) +
  geom_line(aes(linetype = method), 
  lwd = 0.5, data = rarefac.line.bac) +
  scale_colour_discrete(guide = F) +
  scale_x_continuous(limits = c(0,1e+5)) +
  facet_grid(~treatment)+ 
  labs(x="Sample size", y="Species richness") +
  theme_bw() + 
  theme(legend.position="none",
        axis.ticks = element_blank())

richness.euk <- ggplot(rarefac.euk, 
  aes(x=x, y=y, colour=site)) +
  geom_line(aes(linetype = method), 
  lwd = 0.5, data = rarefac.line.euk)+
  scale_colour_discrete(guide = F) +
  scale_x_continuous(limits = c(0,1e+5)) +
  facet_grid(~treatment)+ 
  labs(x="Sample size", y="Species richness") +
  theme_bw() + 
  theme(legend.position="none",
        axis.ticks = element_blank())

###################################

## CREATE SUMMARIES ##

richness <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Species richness",]
simpson <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Simpson diversity",]

# Calculate Shannon Index & evenness
shannon <- estimate_richness(
  pseq.bac.abs, measures=c("Shannon")) %>%
  rownames_to_column("internal_id") %>%
  arrange(internal_id)

bac.div <- data.frame(
  internal_id = iNEXT.bac$DataInfo$site,
  richness = richness$Observed,
  simpson = simpson$Observed,
  shannon = shannon$Shannon) %>%
  mutate(evenness = shannon/log(richness))  %>%
  tibble::add_column(
    type="BAC")
  

###################################

# Combine EUK data
richness <- iNEXT.euk$AsyEst[
  iNEXT.euk$AsyEst$Diversity=="Species richness",]
simpson <- iNEXT.euk$AsyEst[
  iNEXT.euk$AsyEst$Diversity=="Simpson diversity",]

# Calculate Shannon Index & evenness
shannon <- estimate_richness(
  pseq.euk.abs, measures=c("Shannon")) %>%
  rownames_to_column("internal_id") %>%
  arrange(internal_id)

euk.div <- data.frame(
  internal_id = iNEXT.euk$DataInfo$site,
  richness = richness$Observed,
  simpson = simpson$Observed,
  shannon = shannon$Shannon) %>%
  mutate(evenness = shannon/log(richness)) %>%
  tibble::add_column(
    type="EUK") 

###################################

# Merge data
AlphaDiv <- rbind(bac.div, euk.div) %>%
  left_join(dplyr::select(ENV, 
   internal_id, time, treatment, type, extraction), 
   by=c("internal_id"="internal_id","type"="type")) 

# Plot curves
plot_grid(
  richness.bac, 
  richness.euk, 
  coverage.bac, 
  coverage.euk,  
labels = c("BAC","EUK","BAC","EUK"),
label_x = 0.1,
label_size = 9,
ncol = 2)

###################################

# remove temp-data
rm(richness, simpson, shannon)
rm(list = ls(pattern =
   "cover.point.*|cover.line.*|rarefac.point.*|rarefac.line.*"))
