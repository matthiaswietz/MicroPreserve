
############################################################################################
   ###  Fixative Experiment - amplicon analysis  ###
############################################################################################

# This script: load / format amplicons and metadata

library(gtools)
library(dplyr)
library(phyloseq)
library(ampvis2)


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
  row.names = 1)

euk.tax <- read.table(
  "euk_tax.txt",
  h = T, 
  sep = "\t", 
  stringsAsFactors = F, 
  row.names = 1)

# remove chloroplast (BAC)
# remove Mammalia/Vertebrata/arthropoda (EUK)
bac.tax <- bac.tax[-grep(
  "Chloroplast|unknown_unclassified", bac.tax$Order),]
euk.tax <- euk.tax[-grep(
  "Mammalia|Vertebrata|Arthropoda", euk.tax$Phylum),]

# Reformat to match
bac.otu <- merge(
  bac.otu, bac.tax, by="row.names", all.x=F)
bac.otu <- bac.otu[mixedorder(bac.otu$Row.names),]
rownames(bac.otu) = bac.otu$Row.names
bac.otu <- bac.otu[, -c(1, 53:57)]

euk.otu <- merge(euk.otu, euk.tax, by = "row.names", all.x=F)
euk.otu <- euk.otu[mixedorder(euk.otu$Row.names),]
rownames(euk.otu) = euk.otu$Row.names
euk.otu <- euk.otu[, -c(1, 57:61)]

# confirm natural order of rows (OTU1, OTU2 etc)
bac.otu <- bac.otu[mixedorder(rownames(bac.otu)),]
bac.tax <- bac.tax[mixedorder(rownames(bac.tax)),]

euk.otu <- euk.otu[mixedorder(rownames(euk.otu)),]
euk.tax <- euk.tax[mixedorder(rownames(euk.tax)),]

# shorten/modify names
bac.tax <- bac.tax %>%
  mutate(across(everything(),~gsub("Clade","SAR11_Clade", .))) %>%
  mutate(across(everything(),~gsub("Candidatus","Cand", .))) %>%
  mutate(across(everything(),~gsub("Roseobacter_clade_NAC11-7_lineage","Roseobacter_NAC11-7", .))) %>%
  mutate(across(everything(),~gsub("_marine_group","", .))) %>%
  mutate(across(everything(),~gsub("_terrestrial_group","", .))) %>%
  mutate(across(everything(),~gsub("_CC9902","", .))) %>%
  mutate(across(everything(),~gsub("unclassified","uc", .))) %>%
  mutate(across(everything(),~gsub("cl_|_or|_fa|_ge","_uc", .))) %>%
  mutate(across(everything(),~gsub("_uc_uc","_uc", .))) %>%
  mutate(across(everything(),~gsub("(Marine_group_B)","", ., fixed=T))) %>%
  mutate(across(everything(),~gsub("(SAR406_clade)","SAR406", ., fixed=T)))

euk.tax <- euk.tax %>%
  mutate(across(everything(),~gsub("unclassified","uc", .))) %>%
  mutate(across(everything(),~gsub("cl_|_or|_fa|_ge","_uc", .))) %>%
mutate(across(everything(),~gsub("_uc_uc","_uc", .))) 

####################################

# for Ampvis: add "Kingdom" column 
bac.tax$Kingdom = "BAC"
euk.tax$Kingdom = "EUK"

# for Ampvis: add *Species* column (>> copy *Genus*)
bac.tax$Species <- bac.tax$Genus
euk.tax$Species <- euk.tax$Genus

# sort columns for correct order
bac.tax <- bac.tax[, c(
  "Kingdom","Phylum","Class","Order",
  "Family","Genus","Species")]
euk.tax <- euk.tax[, c(
  "Kingdom","Phylum","Class","Order",
  "Family","Genus","Species")]


############################################################################################
  ### METADATA ###
############################################################################################

ENV <- read.table(
  "metadata.txt", 
  h=T, sep="\t", stringsAsFactors=F)

# Generate factors for parameter of interest
ENV$type <- factor(ENV$type, levels = c(
  "BAC","EUK"))
ENV$time <- factor(ENV$time, levels = c(
  "T0","10w","28w","50w"))
ENV$treatment <- factor(ENV$treatment, levels = c(
  "none_ref","HgCl","Formalin"))
ENV$extraction <- factor(ENV$extraction, levels = c(
  "PW","NS"))

# Subset for BAC / EUK
bac.env <- subset(ENV, type=="BAC")
row.names(bac.env) = bac.env$internal_id
euk.env <- subset(ENV, type=="EUK")
row.names(euk.env) = euk.env$internal_id

# Check order so ENV and OTU tables match
bac.env <- bac.env[mixedorder(bac.env$internal_id),]
euk.env <- euk.env[mixedorder(euk.env$internal_id),]


############################################################################################
  ###  PHYLOSEQ  ###
############################################################################################

asv = otu_table(bac.otu, taxa_are_rows=T)
tax = tax_table(bac.tax)
rownames(tax) <- rownames(asv)

pseq.bac.abs <- phyloseq(
  otu_table(asv, taxa_are_rows=F), 
  sample_data(bac.env), 
  tax_table(tax)) %>%
  filter_taxa(function(x) 
    sum(x>3) > (0.03*length(x)), T)

# Fix tax-IDs
colnames(pseq.bac.abs@tax_table) <- c(
 "Kingdom","Phylum","Class","Order",
 "Family","Genus","Species")

####################################

asv = otu_table(euk.otu, taxa_are_rows=T)
tax = tax_table(euk.tax)
rownames(tax) <- rownames(asv)

pseq.euk.abs <- phyloseq(
  otu_table(asv, taxa_are_rows=F), 
  sample_data(euk.env), 
  tax_table(tax)) %>%
  filter_taxa(function(x) 
    sum(x>3) > (0.03*length(x)), T)

# Fix tax-IDs
colnames(pseq.euk.abs@tax_table) <- c(
  "Kingdom","Phylum","Class","Order",
  "Family","Genus","Species")

####################################

## transform to rel. abundance

pseq.bac.rel = transform_sample_counts(
  pseq.bac.abs, function(x) x / sum(x) * 100) 
pseq.euk.rel = transform_sample_counts(
  pseq.euk.abs, function(x) x / sum(x) * 100) 


############################################################################################
  ###  AMPVIS ###
############################################################################################

ampvis.bac <- data.frame(
  OTU = rownames(
    phyloseq::otu_table(pseq.bac.abs)@.Data),
  phyloseq::otu_table(pseq.bac.abs)@.Data,
  phyloseq::tax_table(pseq.bac.abs)@.Data,
  check.names = F)
bac.env <- data.frame(
  phyloseq::sample_data(pseq.bac.abs), 
  check.names = F) %>%
  mutate(Sample=rownames(.), .before = sample_title)
ampvis.bac <- amp_load(ampvis.bac, bac.env)

ampvis.euk <- data.frame(
  OTU = rownames(
    phyloseq::otu_table(pseq.euk.abs)@.Data),
  phyloseq::otu_table(pseq.euk.abs)@.Data,
  phyloseq::tax_table(pseq.euk.abs)@.Data,
  check.names = F)
euk.env <- data.frame(
  phyloseq::sample_data(pseq.euk.abs), 
  check.names = F) %>%
  mutate(Sample=rownames(.), .before = sample_title)
ampvis.euk <- amp_load(ampvis.euk, euk.env)


############################################################################################
   ### Set plotting colors ###
############################################################################################
  
colors <- c(
  "none_ref" = "#50516b",
  "HgCl" = "#2a7fff", 
  "Formalin" = "#f7bd10")


######################################################################

# remove temporary datasets
rm(asv, tax, euk.env, bac.env)
