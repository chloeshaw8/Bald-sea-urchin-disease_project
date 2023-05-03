# ============================================================
'R code for Microbiome analysis

Christina Pavloudi
cpavloudi@gwu.edu
https://cpavloud.github.io/mysite/

	Copyright (C) 2023 Christina Pavloudi
  
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
  
    This script is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.'

# =============================================================

############################LOAD LIBRARIES #######################################

# List of packages needed
.packages = c("vegan", "ecodist", "GGally", "BiocManager", "ggplot2", "tidyverse", "RColorBrewer")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

packageVersion("vegan")
packageVersion("ecodist")
packageVersion("GGally")
packageVersion("ggplot2")
packageVersion("tidyverse")
packageVersion("RColorBrewer")

BiocManager::install("phyloseq")
install.packages("phyloseq")
library(phyloseq); packageVersion("phyloseq")

# Define a default theme for ggplot graphics
theme_set(theme_bw()) 

########################################################################################
############################# 16S rRNA amplicons #######################################
########################################################################################

############################ Preparation of Data #######################################
#import the OTU table (or else biotic data)
bac <- read.csv("seq_table_noblank.csv", sep = ",", header=TRUE, row.names = 1)

#import the taxonomy table
taxonomybac <- read.csv("seq_Taxonomy_silva_merged.csv", sep = ",", header=FALSE, row.names = 1, na.strings=c("","NA"))

colnames(taxonomybac) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#check where there are NA values in the taxonomy data frame
colSums(is.na(taxonomybac))

taxonomy <- taxonomybac

#fill in Phylum column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Phylum[i])){
      taxonomy$Phylum[i]=taxonomy$Kingdom[i]
    }
  if (sum(is.na(taxonomy$Phylum))==0) {
    break
  }
}

#fill in Class column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Class[i])){
      taxonomy$Class[i]=taxonomy$Phylum[i]
    }
  if (sum(is.na(taxonomy$Class))==0) {
    break
  }
}

#fill in Order column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Order[i])){
      taxonomy$Order[i]=taxonomy$Class[i]
    }
  if (sum(is.na(taxonomy$Order))==0) {
    break
  }
}

#fill in Family column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Family[i])){
      taxonomy$Family[i]=taxonomy$Order[i]
    }
  if (sum(is.na(taxonomy$Family))==0) {
    break
  }
}

#fill in Genus column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Genus[i])){
      taxonomy$Genus[i]=taxonomy$Family[i]
    }
  if (sum(is.na(taxonomy$Genus))==0) {
    break
  }
}

#fill in Species column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Species[i])){
      taxonomy$Species[i]=taxonomy$Genus[i]
    }
  if (sum(is.na(taxonomy$Species))==0) {
    break
  }
}


taxonomybac <- taxonomy

#check where there are NA values in the taxonomy data frame and in the biotic data
colSums(is.na(taxonomybac))
colSums(is.na(bac))

#convert the taxonomy data from data frame to matrix
taxonomy_matrix_bac <- as.matrix(taxonomybac)

# prepare the object for the phyloseq object
TAX_BAC = tax_table(taxonomy_matrix_bac)

#convert the biotic data from data frame to matrix
biotic_matrix_bac <- as.matrix(bac)

#tranpose biotic data for the calculation of diversity indices
biotic_presabs <- decostand(bac, method = "pa")
biotic_presabs_matrix <- as.matrix(biotic_presabs)
OTU_PR = otu_table(biotic_presabs_matrix, taxa_are_rows = TRUE)

#prepare the object for the phyloseq object
OTU_BAC = otu_table(biotic_matrix_bac, taxa_are_rows = TRUE)

#import the metadata of the samples
metadata_physeq <- read.csv("metadata.csv", header=TRUE, row.names = 1) 
# prepare the objects for the phyloseq object
META = sample_data(metadata_physeq)
head(META)

######################## PHYLOSEQ analysis #######################################

# combine them all to create the phyloseq object
physeq_bac = phyloseq(OTU_BAC, TAX_BAC, META)
physeq_PR = phyloseq(OTU_PR, TAX_BAC, META)
physeq = physeq_bac

#Check if there are ASVs with no counts and how many there are
any(taxa_sums(physeq) == 0)
sum(taxa_sums(physeq) == 0)
#Remove ASVs/OTUs/taxa that are empty, ie. that have no counts
physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
# Remove samples that are now empty, ie. that have no counts
physeq <- prune_samples(sample_sums(physeq) > 0, physeq)

#Check if there are ASVs with no counts and how many there are
any(taxa_sums(physeq_PR) == 0)
sum(taxa_sums(physeq_PR) == 0)
#Remove ASVs/OTUs/taxa that are empty, ie. that have no counts
physeq_PR <- prune_taxa(taxa_sums(physeq_PR) > 0, physeq_PR)
# Remove samples that are now empty, ie. that have no counts
physeq_PR <- prune_samples(sample_sums(physeq_PR) > 0, physeq_PR)

#get the data frame from the phyloseq object
pd <- psmelt(physeq)
pd_PR <- psmelt(physeq_PR)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd[,c("Phylum")])))
HowManyPhylaPR <- length(levels(as.factor(pd_PR$Phylum)))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)
PhylaPalettePR = getPalette(HowManyPhylaPR)

pd$Status <- factor(pd$Status,      # Reordering group factor levels
                       levels = c("Healthy", "Diseased", "Recovered"))
pd_PR$Status <- factor(pd_PR$Status,      # Reordering group factor levels
                          levels = c("Healthy", "Diseased", "Recovered"))

#merge the OTUs at the Phylum level
physeq_merged_Phylum <- tax_glom(physeq, "Phylum")
ps0 <- transform_sample_counts(physeq_merged_Phylum, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#save the final phyto otu table with the taxonomies
write.csv(OTU_TAX_merged, "16S_OTU_TAX_Merged_Phylum.csv")

#merge the OTUs at the Class level
physeq_merged_Class <- tax_glom(physeq, "Class")
ps0 <- transform_sample_counts(physeq_merged_Class, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#save the final phyto otu table with the taxonomies
write.csv(OTU_TAX_merged, "16S_OTU_TAX_Merged_Class.csv")

#merge the OTUs at the Species level
physeq_merged_Species <- tax_glom(physeq, "Species")
ps0 <- transform_sample_counts(physeq_merged_Species, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#save the final phyto otu table with the taxonomies
write.csv(OTU_TAX_merged, "16S_OTU_TAX_Merged_Species.csv")

# create the nmds plot colour coded by e.g. Type (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq, ord.nmds.bray, color="Status", shape = "Sample_type") +
  geom_point(size=5)+
  scale_color_manual(values= c("#ae2012","#55a630", "#0077b6"))+
  stat_ellipse()

#p1 <- plot_ordination(physeq, ord.nmds.bray, color="Type", title="Bray NMDS", shape = "Sampling_month", label = "Name")
ggsave("16S_P1_Silva_Status.png", width = 8, height = 8, dpi = 600)

#PERMANOVA in phyloseq
metadata_permanova <- as(sample_data(physeq), "data.frame")
permanova.Status <- adonis2(distance(physeq, method="bray") ~ Status, data = metadata_permanova)
permanova.Status

# plot the diversity indices with colour coding by e.g. dpw (info included in the metadata) 
alpha_meas = c("Observed", "Chao1", "ACE")
p_alpha_meas <- plot_richness(physeq, "Status", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas+
  geom_boxplot(data=p_alpha_meas$data, aes(x=Status, y=value, fill = Status), fatten=NULL)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), 
               width= 0.75, linetype= "solid")+
  scale_fill_manual(values= c("#ae2012", "#55a630", "#0077b6"))+
  geom_jitter(width=0, size= 2, color= "black", aes(shape = Sample_type))+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  labs(title= "Phylogenetic Diversity", y= "Diversity Index", x="Group")+
  scale_x_discrete(limits= c("Diseased", "Recovered", "Healthy"))

res.aov <- aov(value ~ Status, data= p_alpha_meas$data)  
summary(res.aov)



######################### LefSe etc. ###########################################

#For the LefSe
BiocManager::install("microbiomeMarker")
install.packages("microbiomeMarker")
library(microbiomeMarker); packageVersion("microbiomeMarker"); citation("microbiomeMarker")
install.packages("MicrobiotaProcess")
BiocManager::install("MicrobiotaProcess")
library(MicrobiotaProcess); packageVersion("MicrobiotaProcess"); citation("MicrobiotaProcess")
library(patchwork)
#For the upset plot
install.packages("UpSetR")
library(UpSetR); packageVersion("UpSetR"); citation("UpSetR")
install.packages("ComplexUpset")
library(ComplexUpset); packageVersion("ComplexUpset"); citation("ComplexUpset")
# for the kruskal_test and wilcox_test
install.packages("coin")
library(coin)

# three groups
#mg_anova <- run_test_multiple_groups(
#  physeq,
#  group = "Status",
#  method = "anova"
#)
#mg_anova

#LefSe on ASV level
OTU_lefse <- run_lefse(
  physeq,
  group = "Status",
  subgroup = "Sample_type",
  taxa_rank = "none", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
  transform = c("identity"), #return the original data without any transformation (default).
  norm = "CPM", # do not normalize.
  kw_cutoff = 0.05,
  lda_cutoff = 3.1,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = "0",
  sample_min = 7,
  only_same_subgrp = FALSE,
  curv = FALSE
)


markertable_status_type <- marker_table(OTU_lefse)
markertable_status_type <- as.data.frame(markertable_status_type)

OTU_lefse
head(marker_table(OTU_lefse))
write.csv(marker_table(OTU_lefse), "lefse_status_type.csv")


heatmap <- plot_heatmap(OTU_lefse, transform = "log10p", group = "Status", 
                        column_order = c("D1", "D2", "D3", "D4", "SW_D1", "SW_D2", "R1", "R2", "R3", "R4", "SW_R1", "SW_R2", "H1", "H2", "H3", "H4", "SW_H1", "SW_H2"),
                        sample_label = TRUE, 
                        annotation_col = c("Healthy" = "#55a630", "Recovered" = "#0077b6", "Diseased" = "#ae2012"))
ggsave("lefse_heatmap_status_type.png", width = 8, height = 8, dpi = 600)



####################### UPSET PLOT #############################################

### Create upset diagram
upsetda <- get_upset(obj=physeq, factorNames="Status")
#upset(upsetda, sets=unique(as.vector(sample_data(physeq)$Status)), 
#      sets.bar.color = "#56B4E9",
#      order.by = "freq", 
#      empty.intersections = "on",mainbar.y.label	= "Number of ASVs", 
#      point.size = 8, line.size	= 2)
#ggsave("upset_status.png", width = 12, height = 8, dpi = 600)
write.csv(upsetda, "upset.csv")


Status = c("Healthy", "Diseased", "Recovered")

upset(
  upsetda,
  Status,
  name='Status',
  queries=list(
    upset_query(
      intersect=c('Healthy', 'Recovered'),
      color='black',
      fill='black',
      only_components=c('intersections_matrix', 'Intersection size')
    )),
  base_annotations=list(
    'Number of ASVs'=intersection_size(
      text=list(
        vjust=-0.3,
        hjust=--0.4,
        angle=0, size=5
      )
    )
    + annotate(
      geom='text', size=5, x=Inf, y=Inf,
      label=paste('Total ASVs:', nrow(upsetda)),
      vjust=1, hjust=1
    )
  ),
  min_size=5,
  width_ratio=0.1,  stripes='white',
  set_sizes=FALSE,
  themes=upset_default_themes(text=element_text(size=15))
)
ggsave("upset_status.png", width = 12, height = 8, dpi = 600)

##################### Rarefaction #############################################

# for reproducibly random number
#set.seed(1024)
#rareres <- get_rarecurve(obj=physeq, chunks=400)

#p_rare <- ggrarecurve(obj=rareres,
#                      indexNames=c("Observe","Chao1","ACE"),
#) +
#  theme(legend.spacing.y=unit(0.01,"cm"),
#        legend.text=element_text(size=4))

#prare1 <- ggrarecurve(obj=rareres, factorNames="Status",
#                      indexNames=c("Observe", "Chao1", "ACE")
#) +
#  scale_fill_manual(values = c("Healthy" = "green", "Recovered" = "blue",
#                               "Diseased" = "pink"))+
#  scale_color_manual(values = c("Healthy" = "green", "Recovered" = "blue",
#                                "Diseased" = "pink"))+
#  theme_bw()+
#  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
#        strip.background = element_rect(colour=NA,fill="grey"),
#        strip.text.x = element_text(face="bold"))          
#
#prare2 <- ggrarecurve(obj=rareres,
#                      factorNames="Status",
#                      shadow=FALSE,
#                     indexNames=c("Observe", "Chao1", "ACE")
#) +
#  scale_color_manual(values = c("Healthy" = "green", "Recovered" = "blue",
#                                "Diseased" = "pink"))+
#  theme_bw()+
#  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
#        strip.background = element_rect(colour=NA,fill="grey"),
#        strip.text.x = element_text(face="bold"))

#p_rare / prare1 / prare2

MPSE %<>%
  mp_rrarefy(.abundance=Abundance) %>%
  mp_cal_rarecurve(.abundance=RareAbundance, chunks=500)
p_rare <- MPSE %>%
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve,
    .alpha = c(Observe, Chao1, ACE),
  ) +
  theme(
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.01,"cm"),
    legend.text = element_text(size=4)
  )
prare1 <- MPSE %>%
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve,
    .alpha = c(Observe, Chao1, ACE),
    .group = Status
  ) +
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  theme_bw()+
  theme(
    axis.text=element_text(size=8), panel.grid=element_blank(),
    strip.background = element_rect(colour=NA,fill="grey"),
    strip.text.x = element_text(face="bold")
  )
prare2 <- MPSE %>%
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve,
    .alpha = c(Observe, Chao1, ACE),
    .group = Status,
    plot.group = TRUE
  ) +
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols) +
  theme_bw()+
  theme(
    axis.text=element_text(size=8), panel.grid=element_blank(),
    strip.background = element_rect(colour=NA,fill="grey"),
    strip.text.x = element_text(face="bold")
  )
(p_rare / prare1 / prare2) + patchwork::plot_annotation(tag_levels="A")


ggsave("rarefaction.png", width = 12, height = 14, dpi = 600)
