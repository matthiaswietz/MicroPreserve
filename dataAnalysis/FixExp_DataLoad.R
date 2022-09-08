
############################################################################################
   ###  Fixative Experiment - Amplicons  ###
############################################################################################

# This script: load / format amplicons and metadata
# Five reagents are included in the original set of fastqs
# Here we only analyze HgCl, formalin, RNAlater and DNAgard
# PBuffer is not evaluated (only a stabilizing PO4 buffer)

setwd("/AWI_MPI/FRAM/FixativeExp/Rstats")
load("FixExp.Rdata")

# Load packages and colors
library(gtools)
library(dplyr)
library(phyloseq)
library(ampvis2)
library(ggplot2)


############################################################################################
   ### OTU IMPORT ###
############################################################################################

bac.otu <- read.table(
  "bac_otu.txt",
  h = T,
  sep = "\t",
  check.names = F, 
  row.names = 1)

euk.otu <- read.table(
  "euk_otu.txt",
  h = T,
  sep = "\t",
  check.names = F, 
  row.names = 1)


############################################################################################
   ### TAXONOMY ###
############################################################################################

bac.tax <- read.table(
  "bac_tax.txt",
  h = T, 
  sep = "\t", 
  stringsAsFactors = F, 
  row.names = 1) %>%
  na_if("")

euk.tax <- read.table(
  "euk_tax.txt",
  h = T, 
  sep = "\t", 
  stringsAsFactors = F, 
  row.names = 1) %>%
  na_if("")

# remove chloroplast (BAC)
# remove Metazoa (EUK)
bac.tax <- bac.tax[-grep(
  "Chloroplast|unknown_unclassified", bac.tax$Order),]
euk.tax <- euk.tax[-grep(
  "Metazoa", euk.tax$Phylum),]

# Reformat to match
bac.otu <- merge(
  bac.otu, bac.tax, by="row.names", all.x=F)
bac.otu <- bac.otu[mixedorder(bac.otu$Row.names),]
rownames(bac.otu) = bac.otu$Row.names
bac.otu <- bac.otu[, -c(1, 142:147)]

euk.otu <- merge(euk.otu, euk.tax, by = "row.names", all.x=F)
euk.otu <- euk.otu[mixedorder(euk.otu$Row.names),]
rownames(euk.otu) = euk.otu$Row.names
euk.otu <- euk.otu[, -c(1,140:146)]

# Confirm natural row order (OTU1, OTU2 etc)
bac.otu <- bac.otu[mixedorder(rownames(bac.otu)),]
bac.tax <- bac.tax[mixedorder(rownames(bac.tax)),]

euk.otu <- euk.otu[mixedorder(rownames(euk.otu)),]
euk.tax <- euk.tax[mixedorder(rownames(euk.tax)),]

######################################

# Add *Species* dummy column (for ampvis)
bac.tax$Species <- bac.tax$Genus

# shorten/modify names
bac.tax <- bac.tax %>%
  mutate(across(everything(),~gsub("Clade","SAR11_Clade", .))) %>%
  mutate(across(everything(),~gsub("Candidatus","Cand", .))) %>%
  mutate(across(everything(),~gsub("Roseobacter_clade_NAC11-7_lineage","Roseobacter_NAC11-7", .))) %>%
  mutate(across(everything(),~gsub("_marine_group","", .))) %>%
  mutate(across(everything(),~gsub("_terrestrial_group","", .))) %>%
  mutate(across(everything(),~gsub("_CC9902","", .))) %>%
  mutate(across(everything(),~gsub("unclassified","uc", .))) %>%
  mutate(across(everything(),~gsub("uncultured", NA, .))) %>%
 # mutate(across(everything(),~gsub("cl_|_or|_fa|_ge","_uc", .))) %>%
  #mutate(across(everything(),~gsub("_uc_uc","_uc", .))) %>%
  mutate(across(everything(),~gsub("(Marine_group_B)","", ., fixed=T))) %>%
  mutate(across(everything(),~gsub("(SAR406_clade)","SAR406", ., fixed=T)))

# Rename BAC-NAs with last known taxrank + "uc"
k <- ncol(bac.tax)-1
for (i in 2:k) {
  if (sum(is.na(bac.tax[, i])) >1) {
    temp <- bac.tax[is.na(bac.tax[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], " uc", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    bac.tax[is.na(bac.tax[, i]), ] <- temp}
  if (sum(is.na(bac.tax[, i]))==1) {
    temp <- bac.tax[is.na(bac.tax[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], " uc", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    bac.tax[is.na(bac.tax[, i]),] <- temp
  }
}
bac.tax[is.na(bac.tax[, (k+1)]), (k+1)] <- paste(
  bac.tax[is.na(bac.tax[, (k+1)]), k], " uc", sep="")

######################################

# Rename EUK-NAs with last known taxrank + "uc"
k <- ncol(euk.tax)-1
for (i in 2:k) {
  if (sum(is.na(euk.tax[, i])) >1) {
    temp <- euk.tax[is.na(euk.tax[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], " uc", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    euk.tax[is.na(euk.tax[, i]), ] <- temp}
  if (sum(is.na(euk.tax[, i]))==1) {
    temp <- euk.tax[is.na(euk.tax[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], " uc", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    euk.tax[is.na(euk.tax[, i]),] <- temp
  }
}
euk.tax[is.na(euk.tax[, (k+1)]), (k+1)] <- paste(
  euk.tax[is.na(euk.tax[, (k+1)]), k], " uc", sep="")
                             
euk.tax <- euk.tax %>%
  mutate(across(everything(),~gsub("_unclassified"," uc", .))) %>%
  mutate(across(everything(),~gsub("Dino-Group-I-Clade-","Dino-I-", .))) %>%
  mutate(across(everything(),~gsub("Dino-Group-II-Clade-","Dino-II-", .))) %>%
  mutate(across(everything(),~gsub("Polar-centric-","", .))) %>%
  mutate(across(everything(),~gsub("Radial-centric-basal-","", .))) %>%
  mutate(across(everything(),~gsub("Chrysophyceae_Clade-","Chrysophyceae ", .))) %>%
  mutate(across(everything(),~gsub("Stephanoecidae_Group_","Stephanoecidae ", .))) %>%
  mutate(across(everything(),~gsub("Pirsonia_Clade_","Pirsonia ", .))) %>%
  mutate(across(everything(),~gsub("_X|_XX|_XXX|_XXXX"," uc", .)))  %>%
  #mutate(across(everything(),~gsub("cl_|_or|_fa|_ge","_uc", .))) %>%
 mutate(across(everything(),~gsub(" uc uc"," uc", .))) 


############################################################################################
  ### METADATA ###
############################################################################################

ENV <- read.table(
  "metadata.txt", 
  h=T, sep="\t", 
  stringsAsFactors=F) %>%
  na_if("")

# Generate factors for parameter of interest
ENV$locus_tag <- factor(ENV$locus_tag, levels = c(
  "16S","18S"))
ENV$time <- factor(ENV$time, levels = c(
  "T0","10w","28w","50w"))
ENV$preservation <- factor(ENV$preservation, levels = c(
  "none_ref","HgCl","Formalin",
  "RNAlater","DNAgard","PBuffer"))
ENV$extraction <- factor(ENV$extraction, levels = c(
  "PowerWater","NucleoSpin"))

#################################

# Import Gfbio sample_titles
names.bac <- read.table("names_bac.txt", header=T)
names.euk <- read.table("names_euk.txt", header=T)

# Sort naturally
bac.otu <- bac.otu[,mixedsort(names(bac.otu))]
euk.otu <- euk.otu[,mixedsort(names(euk.otu))]
names.bac <- names.bac[mixedorder(names.bac$swarm_id),]
names.euk <- names.euk[mixedorder(names.euk$swarm_id),]

# Rename
names(bac.otu) = names.bac$internal_id
names(euk.otu) <- names.euk$internal_id

#################################

# Subset for BAC / EUK
bac.env <- subset(ENV, locus_tag=="16S")
row.names(bac.env) = bac.env$internal_id
euk.env <- subset(ENV, locus_tag=="18S")
row.names(euk.env) = euk.env$internal_id

# Confirm consistent order 
bac.otu <- bac.otu[,mixedsort(names(bac.otu))]
euk.otu <- euk.otu[,mixedsort(names(euk.otu))]
bac.env <- bac.env[mixedorder(bac.env$internal_id),]
euk.env <- euk.env[mixedorder(euk.env$internal_id),]

# Exclude failed samples from metadata
bac.env <- bac.env[
  bac.env$internal_id %in% names(bac.otu), ] 
euk.env <- euk.env[
  euk.env$internal_id %in% names(euk.otu), ] 


############################################################################################
  ###  PHYLOSEQ  ###
############################################################################################

asv = otu_table(bac.otu, taxa_are_rows=T)
tax = tax_table(bac.tax)
rownames(tax) <- rownames(asv)

bac.abs <- phyloseq(
  otu_table(asv, taxa_are_rows=F), 
  sample_data(bac.env), 
  tax_table(tax)) %>%
  filter_taxa(function(x) 
    sum(x>3) > (0.03*length(x)), T)

# Fix tax-IDs
colnames(bac.abs@tax_table) <- c(
 "Kingdom","Phylum","Class","Order",
 "Family","Genus","Species")

####################################

asv = otu_table(euk.otu, taxa_are_rows=T)
tax = tax_table(euk.tax)
rownames(tax) <- rownames(asv)

euk.abs <- phyloseq(
  otu_table(asv, taxa_are_rows=F), 
  sample_data(euk.env), 
  tax_table(tax)) %>%
  filter_taxa(function(x) 
    sum(x>3) > (0.03*length(x)), T)

# Fix tax-IDs
colnames(euk.abs@tax_table) <- c(
  "Kingdom","Phylum","Class","Order",
  "Family","Genus","Species")

###################################

## transform to rel. abundance

bac.rel = transform_sample_counts(
  bac.abs, function(x) x / sum(x) * 100) 
euk.rel = transform_sample_counts(
  euk.abs, function(x) x / sum(x) * 100) 


############################################################################################
  ###  AMPVIS ###
############################################################################################

ampvis.bac <- data.frame(
  OTU = rownames(
    phyloseq::otu_table(bac.abs)@.Data),
  phyloseq::otu_table(bac.abs)@.Data,
  phyloseq::tax_table(bac.abs)@.Data,
  check.names = F)
env <- data.frame(
  phyloseq::sample_data(bac.abs), 
  check.names = F) %>%
  mutate(Sample=rownames(.), .before=sample_title)
ampvis.bac <- amp_load(ampvis.bac, env)

ampvis.euk <- data.frame(
  OTU = rownames(
    phyloseq::otu_table(euk.abs)@.Data),
  phyloseq::otu_table(euk.abs)@.Data,
  phyloseq::tax_table(euk.abs)@.Data,
  check.names = F)
env <- data.frame(
  phyloseq::sample_data(euk.abs), 
  check.names = F) %>%
  mutate(Sample=rownames(.), .before=sample_title)
ampvis.euk <- amp_load(ampvis.euk, env)


############################################################################################
   ### Set plotting options ###
############################################################################################
  
colors.all <- c(
  "none_ref" = "#50516b",
  "HgCl" = "#2a7fff", 
  "Formalin" = "#f7bd10",
  "RNAlater" = "seagreen4",
  "DNAgard" = "navajowhite3")

colors <- c(
  "none_ref" = "#50516b",
  "HgCl" = "#2a7fff", 
  "Formalin" = "#f7bd10")

theme_classic2 <- function(
  base_size=12, base_family=""){
  theme_bw(base_size=base_size, base_family=base_family) %+replace%
    theme(panel.border = element_blank(),
          axis.line = element_line(colour="black"),
          panel.grid.major = element_line(),
          #panel.grid.major.x = element_blank(),
          #panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(colour="black", size = 0.5),
          legend.key = element_blank()
    )
}


######################################################################

# remove temporary datasets
rm(asv, tax, env, i, j, k)


