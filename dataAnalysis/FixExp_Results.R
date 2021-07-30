
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
   ###  DNA YIELDS  ###
###################################################################################

ggplot(ENV, 
  aes(x=time, y=DNA_ng_µL, fill=treatment)) +
  geom_boxplot() + 
  scale_fill_manual(values=colors) +
  facet_grid(
    extraction~treatment, scales="free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank())


###################################################################################
   ###  ORDINATION  ###
###################################################################################

## BAC 
amp_ordinate(
  ampvis.bac,
  filter_species = 1, 
  type = "PCoA",
  distmeasure = "bray",
  transform = "hellinger", 
  x_axis = 1, 
  y_axis = 2, 
  print_caption = F,
  sample_color_by = "treatment", 
  sample_shape_by = "time",
  #sample_label_by = "sample_title",
  sample_point_size = 5) + 
  scale_shape_manual(values=c(15,18,17,19)) + 
  scale_color_manual(values=c(colors)) + 
  #guides(fill=F, colour=F) +  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.ticks = element_blank())

## EUK
amp_ordinate(
  ampvis.euk,  
 filter_species = 1, 
 type = "PCoA",
 distmeasure = "bray",
 transform = "hellinger", 
 x_axis = 1, 
 y_axis = 2, 
 print_caption = F,
 sample_color_by = "treatment", 
 sample_shape_by = "time",
 sample_point_size = 5) + 
  scale_shape_manual(values=c(15,18,17,19)) + 
  scale_color_manual(values=c(colors)) + 
  #guides(fill=F, colour=F) +  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.ticks = element_blank())


###################################################################################
   ###  AlphaDiversity  ###
###################################################################################

with(subset(AlphaDiv, 
    type==c("BAC") & 
    extraction==c("NS")),
  kwAllPairsDunnTest(
    richness ~ time))

#EUK 
ggplot(
  subset(AlphaDiv, type=="EUK"),
  aes(x=time, y=richness, fill=treatment)) +
  geom_boxplot() + 
  scale_fill_manual(values=colors) +
  facet_grid(extraction~treatment, scales="free") +
  ylab("species richness") +
  ggtitle("Eukaryotes") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank())

ggplot(subset(
  AlphaDiv, type=="EUK"),
aes(x=time, y=simpson, fill=treatment)) +
geom_boxplot() + 
scale_fill_manual(values=colors) +
facet_grid(extraction~treatment, scales = "free") +
ylab("invSimpson index") +
ggtitle("Eukaryotes") +
theme_bw() +
theme(
  panel.grid.minor = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank())

######################

with(subset(
  AlphaDiv, type==c("BAC") & extraction==c("PW")),
  kwAllPairsDunnTest(
    simpson ~ treatment))

ggplot(subset(
  AlphaDiv, type=="BAC"),
aes(x=time, y=richness, fill=treatment)) +
geom_boxplot() + 
scale_fill_manual(values=colors) +
facet_grid(extraction~treatment, scales = "free") +
ylab("species richness") +
ggtitle("Bacteria") +
theme_bw() +
theme(
  axis.ticks = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_blank())

ggplot(subset(
  AlphaDiv, type=="BAC"),
aes(x=time, y=simpson, fill=treatment)) +
geom_boxplot() + 
scale_fill_manual(values=colors) +
facet_grid(extraction~treatment, scales = "free") +
ylab("invSimpson index") +
ggtitle("Bacteria") +
theme_bw() +
theme(
  axis.ticks = element_blank(),
  panel.grid.minor = element_blank(),
  axis.title.x = element_blank())


###################################################################################
   ###  BetaDiv  ###
###################################################################################

envtype1 <- get_variable(pseq.bac.abs, "treatment") 
envtype2 <- get_variable(pseq.euk.abs, "treatment") 

## Map sample type to color 

tipColor1 <- col_factor(
  colors, levels= levels(envtype1))(envtype1) 
tipColor2 <- col_factor(
  colors, levels= levels(envtype2))(envtype2) 

## Hellinger-transform
hel.bac <- transform_sample_counts(
  pseq.bac.rel, function(x) sqrt(x / sum(x)))
hel.euk <- transform_sample_counts(
  pseq.euk.rel, function(x) sqrt(x / sum(x)))

dist.bac <- phyloseq::distance(hel.bac, method="jaccard")
clust.bac <- as.phylo(hclust(dist.bac, method="complete")) 

dist.euk <- phyloseq::distance(hel.euk, method="jaccard")
clust.euk <- as.phylo(hclust(dist.euk, method="complete")) 

# shorten IDs for nicer plotting
clust.bac$tip.label <- gsub('bactV4V5_', '', clust.bac$tip.label) 
clust.bac$tip.label <- gsub('Mercury_Chloride', 'HgCl', clust.bac$tip.label) 
clust.bac$tip.label <- gsub('_weeks', 'w', clust.bac$tip.label) 
clust.bac$tip.label <- gsub('Reference_none', 'T0', clust.bac$tip.label) 

clust.euk$tip.label <- gsub('eukV4V5_', '', clust.euk$tip.label) 
clust.euk$tip.label <- gsub('Mercury_Chloride', 'HgCl', clust.euk$tip.label) 
clust.euk$tip.label <- gsub('_weeks', 'w', clust.euk$tip.label) 
clust.euk$tip.label <- gsub('Reference_none', 'T0', clust.euk$tip.label) 

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

plot.tax1 <- amp_heatmap(
  ampvis.bac,
 tax_aggregate = "Class",
 tax_show = c(
   "Gammaproteobacteria","Alphaproteobacteria",
   "Bacteroidia","Bacteria_uc",
   "Planctomycetacia","Deltaproteobacteria",
   "Actinobacteria","Bacilli"),
 normalise = T,
 plot_values = T,
 round = 0,
 max_abundance = 50,
 min_abundance = 1,
 plot_colorscale = "sqrt",
 plot_legendbreaks = c(5,10,30,50),
 facet_by = c("treatment","extraction"),
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
theme(axis.text = element_text(size=9),
      axis.ticks = element_blank())

plot.tax2 <- amp_heatmap(ampvis.euk, 
  tax_aggregate = "Class",
  tax_show = 11,
  normalise = T,
  plot_values = F,
  round = 0,
  max_abundance = 50,
  min_abundance = 1,
  plot_colorscale = "sqrt",
  plot_legendbreaks = c(5,10,30,50),
  facet_by = c("treatment","extraction"),
  group_by = "time",
  color_vector = c(values =rev(
    scico::scico(4, palette='tokyo')))) +
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
theme(axis.text = element_text(size=9),
    axis.ticks = element_blank())

plot_grid(
  plot.tax1,
  plot.tax2,
  ncol=1,
  align="h",
  axis="lr")

#########################################################

## SUPPLEMENT -- GENUS LEVEL

amp_heatmap(ampvis.bac,
  tax_aggregate = "Genus",
  tax_show = 15,
  normalise = T,
  plot_values = T,
  round = 0,
  max_abundance = 35,
  min_abundance = 1,
  plot_colorscale = "sqrt",
  plot_legendbreaks = c(5,10,20,50),
  facet_by = c("treatment","extraction"),
  group_by = "time",
  color_vector = c(values =rev(
    scico(4, palette='tokyo')))) +
geom_text(aes(
  label = round(Abundance),
  color = ifelse(
    Abundance < 9.5,
    "black","white"))) +
scale_color_identity()

amp_heatmap(
  ampvis.euk,
  tax_aggregate = "Genus",
  tax_show = 15,
  normalise = T,
  plot_values = T,
  round = 0,
  max_abundance = 35,
  min_abundance = 1,
  plot_colorscale = "sqrt",
  plot_legendbreaks = c(5,10,20,50),
 facet_by = c("treatment","extraction"),
 group_by = "time",
 color_vector = c(values =rev(
   scico(4, palette='tokyo')))) +
geom_text(aes(
  label = round(Abundance),
  color = ifelse(
    Abundance < 9.5,
    "black","white"))) +
scale_color_identity() 

###############################

## ALL REPLICATES BAC - no outlier!
amp_heatmap(
  ampvis.bac, 
  tax_aggregate = "Class",
  tax_show = c(
    "Gammaproteobacteria",
    "Alphaproteobacteria",
    "Bacteroidia","Bacteria_uc",
    "Planctomycetacia",
    "Deltaproteobacteria",
    "Actinobacteria","Bacilli"),
  normalise = T,
  plot_values = F,
  round = 0,
  max_abundance = 35,
  min_abundance = 1,
  plot_colorscale = "sqrt",
  plot_legendbreaks = c(5,10,20,50),
  facet_by = c("treatment","extraction"),
  group_by = "sample_title",
  color_vector = c(values =rev(
    scico(4, palette='tokyo')))) 

## ALL REPLICATES EUK - show outlier
amp_heatmap(
  ampvis.euk, 
  tax_aggregate = "Class",
  tax_show = 13,
  normalise = T,
  plot_values = F,
  round = 0,
  max_abundance = 35,
  min_abundance = 1,
  plot_colorscale = "sqrt",
  plot_legendbreaks = c(5,10,20,50),
  facet_by = c("treatment","extraction"),
  group_by = "sample_title",
  color_vector = c(values =rev(
    scico(4, palette='tokyo')))) 


############################################################################################
   ###  Hellinger OTUs  ###
############################################################################################

bac <- pseq.bac.rel %>%
  tax_glom(., "Genus") %>%
  filter_taxa(., function(x) mean(x) > 0.1, T)
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
group_by(variable, treatment, time, 
        extraction, Genus) %>%
summarize_at(c("value"), sum) %>%
ungroup %>%
mutate(variable=gsub(
  'FE[[:digit:]]+_',"", variable)) %>%
group_by(variable, treatment, time, 
        extraction, Genus) %>%
summarize_at(c("value"), mean) %>%
  add_column("type"="bac")

#####################

euk <- pseq.euk.rel %>%
  tax_glom(., "Genus") %>%
  filter_taxa(., function(x) mean(x) > 0.1, T)
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
group_by(variable, treatment, time, 
        extraction, Genus) %>%
summarize_at(c("value"), sum) %>%
ungroup %>%
mutate(variable=gsub(
  'FE[[:digit:]]+_',"", variable)) %>%
group_by(variable, treatment, time, 
        extraction, Genus) %>%
summarize_at(c("value"), mean) %>%
add_column("type"="euk")

###################

ggplot(data=subset(
  hel1, Genus %in% c(
    "Colwellia","SAR11_Clade_Ia",
    "Pseudohongiella","Sulfitobacter",
    "Amphritea","Pseudoalteromonas",
    "Flavobacteriaceae_uc")
  & treatment!="Formalin")) + 
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

ggplot(data=subset(
  hel2, Genus %in% c(
    "Cercozoa_uc","Eukaryota_uc",
    "MAST-2_uc","Picomonas",
    "Stephanoecidae_uc_uc",
    "Thalassiosira","Paraphysomonas",
    "Ochrophyta_uc","Pirsonia") & treatment!="Formalin")) + 
aes(x=time, y=value, fill=Genus) +
geom_bar(stat="identity",position="stack") +
scale_fill_manual(
  values = alpha(scico(9, palette='roma'), 0.7)) +
facet_grid(~extraction) +
scale_y_continuous(
  expand=c(0.02,0.02),
  breaks=c(0,0.5,1)) +
ylab("Hellinger-transformed abundance") +
theme_bw() + 
theme(#axis.text.x = element_text(angle=90),
      axis.title.x=element_blank(),
      axis.ticks = element_blank())


############################################################################################
   ###  NUMBER OF OTUs  ###
############################################################################################

count1 <- subset_taxa(
  pseq.euk.rel, Class=="Diatomea") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "treatment","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Diatomea')

count2 <- subset_taxa(
  pseq.euk.rel, Phylum=="Cercozoa") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "treatment","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Cercozoa')

count3 <- subset_taxa(
  pseq.euk.rel, Class=="Trebouxiophyceae") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "treatment","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Trebouxiophyceae')

count4 <- subset_taxa(
  pseq.euk.rel, Class=="Thecofilosea") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "treatment","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Thecofilosea')

count5 <- subset_taxa(
  pseq.euk.rel, Class=="Chrysophyceae") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "treatment","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Chrysophyceae')

count6 <- subset_taxa(
  pseq.euk.rel, Class=="Acanthoecida") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "treatment","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Acanthoecida')

count7 <- subset_taxa(
  pseq.euk.rel, Class=="Picomonadida") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "treatment","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Picomonadida')

count8 <- subset_taxa(
  pseq.euk.rel, Phylum=="MAST-2") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "treatment","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='MAST-2')

count.euk <- rbind(
  count1,count2,count3,count4,
  count5,count6,count7,count8)

#########################

count9 <- subset_taxa(
  pseq.bac.abs, Class=="Alphaproteobacteria") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "treatment","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Alphaproteobacteria')

count10 <- subset_taxa(
  pseq.bac.rel, Class=="Gammaproteobacteria") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "treatment","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Gammaproteobacteria')

count11 <- subset_taxa(
  pseq.bac.rel, Class=="Bacteroidia") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "treatment","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Bacteroidia')

count.bac <- rbind(
  count9,count10,count11)

#########################

ggplot(
  count.euk, 
aes(x=factor(time), y=freq, 
    group=extraction, color=extraction)) +
geom_point(size=3) +
geom_line(size=1) +
scale_x_discrete(expand = c(0.08,0.08)) +
scale_y_continuous(n.breaks=4) +
ylab("Sequence number") +
facet_grid(
  tax~treatment, scales="free") +
scale_color_manual(values=c(
  "PW"="#09d6c8",
  "NS"="#9e3f1c")) +
theme_bw()+
theme(axis.text.y = element_text(size=10),
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(
        color = "white", size = 0.01),
      axis.ticks = element_blank()) 

ggplot(
  count.bac, 
aes(x=factor(time), y=freq, 
    group=extraction, color=extraction)) +
geom_point(size=3) +
geom_line(size=1) +
scale_x_discrete(expand = c(0.08,0.08)) +
scale_y_continuous(n.breaks=4) +
ylab("Sequence number") +
facet_grid(
  tax~treatment, scales="free") +
scale_color_manual(values=c(
  "PW"="#09d6c8",
  "NS"="#854d22")) +
theme_bw()+
theme(axis.text.y = element_text(size=10),
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(
        color = "white", size = 0.01),
      axis.ticks = element_blank()) 

##############################################

# remove temp-data
rm(count1,count2,count3,count4,
   count5,count6,count7,count8,
   count9,count10,count11)
