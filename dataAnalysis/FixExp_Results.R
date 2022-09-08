
############################################################################################
   ###  Fixative Experiment - amplicon analysis  ###
############################################################################################

# This script: plotting of results 

library(phyloseq)
library(ampvis2)
library(vegan)
library(ggplot2)
library(scico)
library(dplyr)
library(tibble)
library(tidyr)
library(data.table)
library(PMCMRplus)
library(ape)
library(scales)


###################################################################################
   ###  ORDINATION  ###
###################################################################################

## BAC 
amp_ordinate(amp_subset_samples(
    ampvis.bac, preservation %in% c(
    "none_ref","HgCl","Formalin")),
  filter_species = 1, 
  type = "PCoA",
  distmeasure = "bray",
  transform = "hellinger", 
  x_axis = 1, 
  y_axis = 2, 
  print_caption = F,
  sample_color_by = "preservation", 
  sample_shape_by = "time",
  sample_point_size = 5) + 
  scale_shape_manual(values=c(15,18,17,19)) + 
  scale_color_manual(values=c(colors)) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.ticks = element_blank())

## EUK
amp_ordinate(amp_subset_samples(
  ampvis.euk, preservation %in% c(
    "none_ref","HgCl","Formalin")), 
 filter_species = 1, 
 type = "PCoA",
 distmeasure = "bray",
 transform = "hellinger", 
 x_axis = 1, 
 y_axis = 2, 
 print_caption = F,
 sample_color_by = "preservation", 
 sample_shape_by = "time",
 sample_point_size = 5) + 
  scale_shape_manual(values=c(15,18,17,19)) + 
  scale_color_manual(values=c(colors)) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.ticks = element_blank())


###################################################################################
   ###  AlphaDiversity  ###
###################################################################################

# EUK richness
plot.div1 <- ggplot(subset(
  AlphaDiv, locus_tag=="18S" & preservation %in% c(
    "HgCl","Formalin","none_ref")),
  aes(x=time, y=richness, fill=preservation)) +
  geom_boxplot() + 
  scale_fill_manual(values=colors) +
  facet_grid(extraction~preservation, scales="free") +
  ylab("species richness") +
  ggtitle("Eukaryotes") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank())

# EUK simpson
plot.div2 <- ggplot(subset(
  AlphaDiv, locus_tag=="18S" & preservation %in% c(
    "HgCl","Formalin","none_ref")),
aes(x=time, y=simpson, fill=preservation)) +
geom_boxplot() + 
scale_fill_manual(values=colors) +
facet_grid(extraction~preservation, scales = "free") +
ylab("invSimpson index") +
theme_bw() +
theme(
  panel.grid.minor = element_blank(),
  axis.ticks = element_blank(),
  strip.background.x = element_blank(),
  strip.text.x = element_blank(),
  legend.position = "none",
  axis.title.x = element_blank())

# export size 5.5 x 5
plot_grid(
  plot.div1, 
  plot.div2,
  rel_heights = c(1, 0.9),
  axis="lr",
  align="v",
  ncol = 1)

######################

# BAC richness
plot.div3 <- ggplot(subset(
  AlphaDiv, locus_tag=="16S" & preservation %in% c(
    "HgCl","Formalin","none_ref")),
aes(x=time, y=richness, fill=preservation)) +
geom_boxplot() + 
scale_fill_manual(values=colors) +
facet_grid(extraction~preservation, scales = "free") +
ylab("species richness") +
ggtitle("Bacteria") +
theme_bw() +
theme(
  axis.ticks = element_blank(),
  legend.position = "none",
  axis.text.x = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_blank())

# BAC simpson
plot.div4 <- ggplot(subset(
  AlphaDiv, locus_tag=="16S" & preservation %in% c(
    "HgCl","Formalin","none_ref")),
aes(x=time, y=simpson, fill=preservation)) +
geom_boxplot() + 
scale_fill_manual(values=colors) +
facet_grid(extraction~preservation, scales = "free") +
ylab("invSimpson index") +
theme_bw() +
theme(
  axis.ticks = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background.x = element_blank(),
  strip.text.x = element_blank(),
  legend.position = "none",
  axis.title.x = element_blank())

plot_grid(
  plot.div3, 
  plot.div4,
  rel_heights = c(1, 0.9),
  axis="lr",
  align="v",
  ncol = 1)


###################################################################################
   ###  BetaDiv  ###
###################################################################################

# subset samples; Hellinger-transform
subset.bac <- subset_samples(
  bac.rel, preservation %in% c("HgCl","Formalin","none_ref")) %>%
  transform_sample_counts(function(x) sqrt(x / sum(x)))
subset.euk <- subset_samples(
  euk.rel, preservation %in% c("HgCl","Formalin","none_ref")) %>%
  transform_sample_counts(function(x) sqrt(x / sum(x)))

envtype1 <- get_variable(subset.bac, "preservation") 
envtype2 <- get_variable(subset.euk, "preservation") 

## Map sample type to color 

tipColor1 <- col_factor(
  colors, levels= levels(envtype1))(envtype1) 
tipColor2 <- col_factor(
  colors, levels= levels(envtype2))(envtype2) 

# do once with Jaccard, once with Bray
dist.bac <- phyloseq::distance(subset.bac, method="jaccard")
clust.bac <- as.phylo(hclust(dist.bac, method="complete")) 

dist.euk <- phyloseq::distance(subset.euk, method="jaccard")
clust.euk <- as.phylo(hclust(dist.euk, method="complete")) 

# shorten IDs for nicer plotting
clust.bac$tip.label <- gsub('bactV4V5_', '', clust.bac$tip.label) 
clust.bac$tip.label <- gsub('Mercury_Chloride', 'HgCl', clust.bac$tip.label) 
clust.bac$tip.label <- gsub('_weeks', 'w', clust.bac$tip.label) 
clust.bac$tip.label <- gsub('Reference_none', 'T0', clust.bac$tip.label) 
clust.bac$tip.label <- substring(clust.bac$tip.label, 6)

clust.euk$tip.label <- gsub('eukV4V5_', '', clust.euk$tip.label) 
clust.euk$tip.label <- gsub('Mercury_Chloride', 'HgCl', clust.euk$tip.label) 
clust.euk$tip.label <- gsub('_weeks', 'w', clust.euk$tip.label) 
clust.euk$tip.label <- gsub('Reference_none', 'T0', clust.euk$tip.label) 
clust.euk$tip.label <- substring(clust.euk$tip.label, 6)

# export size 10 x 5
par(mar=c(2,1,1,5))
plot(
  clust.bac, 
  tip.color=tipColor1, 
  direction="rightwards", 
  main="BAC - Jaccard - CompleteLinkage") 
axisPhylo()

plot(
  clust.euk, 
  tip.color=tipColor2, 
  direction="rightwards", 
  main="EUK - Jaccard - CompleteLinkage") 
axisPhylo()


############################################################################################
   ###  TOP TAXA  ###
############################################################################################

## BAC + EUK -- CLASS LEVEL
# export size 3.5 x 6.5

amp_heatmap(amp_subset_samples(
  ampvis.bac, preservation %in% c(
    "none_ref","HgCl","Formalin")),
 tax_aggregate = "Class",
 tax_show = c(
   "Gammaproteobacteria","Alphaproteobacteria",
   "Deltaproteobacteria","Bacteroidia",
   "Planctomycetacia","Bacteria_uc",
   "Actinobacteria","Bacilli"),
 normalise = T,
 plot_values = T,
 round = 0,
 max_abundance = 50,
 min_abundance = 1,
 plot_colorscale = "sqrt",
 plot_legendbreaks = c(5,10,30,50),
 facet_by = c("preservation","extraction"),
 group_by = "time",
color_vector = c(values =rev(
  scico(4, palette='tokyo')))) +
geom_text(aes(
  label = round(Abundance),
  color = ifelse(
    Abundance < 9.5,
    "black","white"))) +
scale_color_identity() +
theme(
  axis.text = element_text(size=11),
  axis.text.x = element_text(angle=0),
  legend.position = "bottom",
  axis.ticks = element_blank())

amp_heatmap(amp_subset_samples(
  ampvis.euk, preservation %in% c(
    "none_ref","HgCl","Formalin")),
  tax_aggregate = "Class",
  tax_show = c(
    "Trebouxiophyceae","Bacillariophyta",
    "Cercozoa uc","Filosa-Thecofilosea",
    "Choanoflagellatea","MAST","Picozoa uc",
    "Filosa-Imbricatea","Dinophyceae",
    "Stramenopiles uc","Spirotrichea"),
  normalise = T,
  plot_values = F,
  round = 0,
  max_abundance = 50,
  min_abundance = 1,
  plot_colorscale = "sqrt",
  plot_legendbreaks = c(5,10,30,50),
  facet_by = c("preservation","extraction"),
  group_by = "time",
  color_vector = c(values =rev(
    scico(4, palette='tokyo')))) +
geom_text(aes(
  label = round(Abundance),
  color = ifelse(
    Abundance < 9.5,
    "black","white"))) +
scale_color_identity() +
theme(
  axis.text = element_text(size=11),
  axis.text.x = element_text(angle=0),
  legend.position = "bottom",
  axis.ticks = element_blank())


############################################################################################
   ###  TOP GENERA - BAC  ###
############################################################################################

amp_heatmap(amp_subset_samples(
  ampvis.bac, preservation %in% c(
    "none_ref","HgCl","Formalin")),
  tax_aggregate = "Genus",
  tax_show = 13,
  normalise = T,
  plot_values = F,
  round = 0,
  max_abundance = 35,
  min_abundance = 1,
  plot_colorscale = "sqrt",
  plot_legendbreaks = c(5,10,20,50),
  facet_by = c("preservation","extraction"),
  group_by = "time",
  color_vector = c(values =rev(
    scico(4, palette='tokyo')))) +
  geom_text(aes(
    label = round(Abundance),
    color = ifelse(
      Abundance < 9.5,
      "black","white"))) +
  scale_color_identity()

#############################

## HELLINGER 
bac <- bac.rel %>%
  tax_glom("Genus") %>%
  filter_taxa(function(x) mean(x) > 0.1, T)
otu_table(bac) = otu_table(decostand(
  otu_table(bac), method="hellinger"), 
  taxa_are_rows=T)

hel1 = as(otu_table(bac), "matrix") %>%
  cbind(as(tax_table(bac), "matrix")) %>%
as.data.frame() %>%
reshape2::melt(id.vars=c(
  "Class","Family","Genus")) %>%
left_join(ENV, by=c(
  "variable"="internal_id")) %>%
drop_na() %>%
mutate(across(value, as.numeric)) %>%
group_by(variable, preservation, time, 
        extraction, Genus) %>%
  summarise(value=sum(value)) %>%
ungroup %>%
mutate(variable=gsub(
  'FE[[:digit:]]+_',"", variable)) %>%
group_by(variable, preservation, time, 
        extraction, Genus) %>%
summarize_at(c("value"), mean) 

# export size 3 x 6
ggplot(data=subset(
  hel1, Genus %in% c(
    "Colwellia","SAR11_Clade_Ia",
    "Pseudohongiella","Sulfitobacter",
    "Amphritea","Pseudoalteromonas",
    "Flavobacteriaceae uc")
  & preservation %in% c("HgCl","none_ref"))) + 
  aes(x=time, y=value, fill=Genus) +
  geom_bar(stat="identity",position="stack") +
  facet_grid(~extraction) +
  scale_fill_manual(
    values = alpha(scico(7, palette='roma'), 0.7)) +
  scale_y_continuous(
    expand=c(0.02,0.02),
    breaks=c(0,0.5,1)) +
  ylab("Hellinger-transformed abundance") +
  theme_bw() +
  theme(#axis.text.x = element_text(angle=90),
    axis.title.x=element_blank(),
    axis.ticks = element_blank())


############################################################################################
###  TOP GENERA - EUK  ###
############################################################################################

amp_heatmap(
  amp_subset_samples(
    ampvis.euk, preservation %in% c(
      "none_ref","HgCl","Formalin")),
  tax_aggregate = "Genus",
  tax_show = 13,
  normalise = T,
  plot_values = F,
  round = 0,
  max_abundance = 35,
  min_abundance = 1,
  plot_colorscale = "sqrt",
  plot_legendbreaks = c(5,10,20,50),
  facet_by = c("preservation","extraction"),
  group_by = "time",
  color_vector = c(values =rev(
    scico(4, palette='tokyo')))) +
  geom_text(aes(
    label = round(Abundance),
    color = ifelse(
      Abundance < 9.5,
      "black","white"))) +
  scale_color_identity() 

#############################

## HELLINGER 
euk <- euk.rel %>%
  tax_glom("Genus") %>%
  filter_taxa(function(x) mean(x) > 0.1, T)
otu_table(euk) = otu_table(decostand(
  otu_table(euk), method="hellinger"), 
  taxa_are_rows=T)

hel2 = as(otu_table(euk), "matrix") %>%
  cbind(as(tax_table(euk), "matrix")) %>%
as.data.frame() %>%
reshape2::melt(id.vars=c(
  "Class","Family","Genus")) %>%
left_join(ENV, by=c(
  "variable"="internal_id")) %>%
drop_na() %>%
mutate(across(value, as.numeric)) %>%
group_by(variable, preservation, time, 
        extraction, Genus) %>%
summarise(value=sum(value)) %>%
ungroup %>%
mutate(variable=gsub(
  'FE[[:digit:]]+_',"", variable)) %>%
group_by(variable, preservation, time, 
        extraction, Genus) %>%
summarize_at(c("value"), mean) 


ggplot(data=subset(
  hel2, Genus %in% c(
    "Reckertia","Stramenopiles uc",
    "MAST-1C uc","Picozoa uc",
    "Stephanoecidae D uc","Cercozoa uc",
    "Mediophyceae uc","Protaspa-lineage uc") & 
   preservation %in% c("HgCl","none_ref"))) +
aes(x=time, y=value, fill=Genus) +
geom_bar(stat="identity",position="stack") +
scale_fill_manual(
  values = alpha(scico(8, palette='roma'), 0.7)) +
facet_grid(~extraction) +
scale_y_continuous(
  expand=c(0.02,0.02),
  breaks=c(0,0.5,1)) +
ylab("Hellinger-transformed abundance") +
theme_bw() + 
theme(#axis.text.x = element_text(angle=90),
      axis.title.x=element_blank(),
      axis.ticks = element_blank())

#############
# remove temp-data
rm(bac,euk)

