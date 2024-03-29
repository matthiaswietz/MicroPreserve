###################################################################################
   ###  All-sample comparison ###
###################################################################################

## BAC community
amp_ordinate(amp_subset_samples(
  ampvis.bac, preservation!="PBuffer"),
  filter_species = 1, 
  type = "PCoA",
  distmeasure = "bray",
  transform = "hellinger", 
  x_axis = 1, 
  y_axis = 2, 
  print_caption = F,
  sample_color_by = "preservation", 
  #sample_shape_by = "time",
  #sample_label_by = "sample_title",
  sample_point_size = 5) + 
  scale_shape_manual(values=c(15,18,17,19)) + 
  scale_color_manual(values=c(colors.all)) + 
  #guides(fill=F, colour=F) +  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank())

###########################

## EUK community
amp_ordinate(amp_subset_samples(
  ampvis.euk, preservation!="PBuffer"),
  filter_species = 1, 
  type = "PCoA",
  distmeasure = "bray",
  transform = "hellinger", 
  x_axis = 1, 
  y_axis = 2, 
  print_caption = F,
  sample_color_by = "preservation", 
  #sample_shape_by = "time",
  sample_point_size = 5) + 
  scale_shape_manual(values=c(15,18,17,19)) + 
  scale_color_manual(values=c(colors.all)) + 
  #guides(fill=F, colour=F) +  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank())


###################################################################################
   ###  DNA extraction yield  ###
###################################################################################

# good fixatives
read.table(
  "DNAconc.txt", h=T, sep="\t") %>%
  filter(preservation %in% c(
    "none_ref","HgCl","Formalin")) %>%
  mutate(extraction=factor(
    extraction, levels=c("PowerWater","NucleoSpin"))) %>%
  mutate(preservation=factor(
    preservation, levels=c(
      "none_ref","HgCl","Formalin"))) %>%
  ggplot(., aes(
    x=time, y=DNA_ng_µL, fill=preservation)) +
  geom_boxplot() + 
  scale_fill_manual(values=colors.all) +
  facet_grid(
    extraction~preservation, scales="free") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")



###################################################################################
   ###  ALL REPLICATES ###
###################################################################################

## BAC
## export size 4x11
amp_heatmap(
  amp_subset_samples(
    ampvis.bac, preservation %in% c(
      "none_ref","HgCl","Formalin")),
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
  facet_by = c("preservation","extraction"),
  group_by = "sample_title",
  color_vector = c(values =rev(
    scico(4, palette='tokyo')))) 

###############################

## EUK 
## export size 5 x 12
amp_heatmap(
  amp_subset_samples(
    ampvis.euk, preservation %in% c(
      "none_ref","HgCl","Formalin")),
  tax_aggregate = "Class",
  tax_show = 13,
  normalise = T,
  plot_values = F,
  round = 0,
  max_abundance = 35,
  min_abundance = 1,
  plot_colorscale = "sqrt",
  plot_legendbreaks = c(5,10,20,50),
  facet_by = c("preservation","extraction"),
  group_by = "sample_title",
  color_vector = c(values =rev(
    scico(4, palette='tokyo')))) 


############################################################################################
   ###  NUMBER OF OTUs  ###
############################################################################################

count1 <- subset_taxa(
  pseq.euk.rel, Class=="Bacillariophyta") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "preservation","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Bacillariophyta')

count2 <- subset_taxa(
  pseq.euk.rel, Phylum=="Cercozoa") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "preservation","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Cercozoa')

count3 <- subset_taxa(
  pseq.euk.rel, Class=="Stramenopiles uc") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "preservation","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Stramenopiles uc')

count4 <- subset_taxa(
  pseq.euk.rel, Class=="Filosa-Imbricatea") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "preservation","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Filosa-Imbricatea')

count5 <- subset_taxa(
  pseq.euk.rel, Class=="Picozoa uc") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "preservation","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Picozoa uc')

count6 <- subset_taxa(
  pseq.euk.rel, Class=="Choanoflagellatea") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "preservation","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Choanoflagellatea')

count7 <- subset_taxa(
  pseq.euk.rel, Class=="MAST") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "preservation","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='MAST')

count8 <- subset_taxa(
  pseq.euk.rel, Class=="Filosa-Thecofilosea") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "preservation","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Filosa-Thecofilosea')

count.euk <- rbind(
  count1,count2,count3,count4,
  count5,count6,count7,count8) %>%
  filter(preservation %in% c(
    "HgCl","Formalin","none_ref"))


#########################

count9 <- subset_taxa(
  pseq.bac.abs, Class=="Alphaproteobacteria") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "preservation","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Alphaproteobacteria')

count10 <- subset_taxa(
  pseq.bac.rel, Class=="Gammaproteobacteria") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "preservation","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Gammaproteobacteria')

count11 <- subset_taxa(
  pseq.bac.rel, Class=="Bacteroidia") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  plyr::count(., c(
    "preservation","extraction","time")) %>%
  tidyr::drop_na() %>%
  mutate(tax='Bacteroidia')

count.bac <- rbind(
  count9,count10,count11) %>%
  filter(preservation %in% c(
    "HgCl","Formalin","none_ref"))

#########################

#expot size 7 x 5
ggplot(count.euk, aes(
    x=factor(time), y=freq, 
    group=extraction, color=extraction)) +
  geom_point(size=3) +
  geom_line(size=1) +
  scale_x_discrete(expand = c(0.08,0.08)) +
  scale_y_continuous(
    n.breaks=4,
    expand = c(0.1,0.1)) +
  ylab("Sequence number") +
  facet_grid(
    tax~preservation, scales="free") +
  scale_color_manual(values=c(
    "PowerWater"="#09d6c8",
    "NucleoSpin"="#854d22")) +
  theme_bw()+
  theme(axis.text.y = element_text(size=10),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background.x = element_rect(
          color = "white", size = 0.01),
        strip.background.y = element_blank(),
        axis.ticks = element_blank()) 

ggplot(count.bac, aes(
    x=factor(time), y=freq, 
    group=extraction, color=extraction)) +
  geom_point(size=3) +
  geom_line(size=1) +
  scale_x_discrete(expand = c(0.08,0.08)) +
  scale_y_continuous(
    n.breaks=4,
    expand = c(0.1,0.1)) +
  ylab("Sequence number") +
  facet_grid(
    tax~preservation, scales="free") +
  scale_color_manual(values=c(
    "PowerWater"="#09d6c8",
    "NucleoSpin"="#854d22")) +
  theme_bw()+
  theme(axis.text.y = element_text(size=10),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background.x = element_rect(
          color = "white", size = 0.01),
        strip.background.y = element_blank(),
        axis.ticks = element_blank()) 

##############################################

# remove temp-data
rm(count1,count2,count3,count4,
   count5,count6,count7,count8,
   count9,count10,count11)
