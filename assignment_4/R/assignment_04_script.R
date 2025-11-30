##***************************
## BINF 6410 Software Tools - Assignment 4
## Assessing Phylogenetic Structure and Community Assembly of Great Lakes Fish
##
## Farah Sadoon 1302190
##
## 2025-12-05
##
##***************************

rm(list = ls())

## _Setting up Libraries --------
library(tidyverse)
library(janitor)
library(picante)
library(ape)
library(stringi)
library(httr)
library(Biostrings)
library(phylotools)
library(ggtree)

##_Loading, Exploring and Cleaning Data ----

#### _ Exploring and Cleaning Great Lakes Fish Data --------
# Downloaded from https://www.glerl.noaa.gov/data/waterlife/glwlexplorer.html 
df_fish <- read_tsv("../data/fish_great_lakes_glerl_noaa.tsv")

# Explore great lakes fish species data frame
class(df_fish)
dim(df_fish)
head(df_fish)
tail(df_fish)
str(df_fish)
names(df_fish)

# Check if any species names are missing
table(df_fish$Species)

# Create a subset of the fish data set that includes only the columns of interest
df_fish2 <- df_fish %>%
  dplyr::select("Species", "Superior", "Michigan", "Huron", "Erie", "Ontario") %>%
  mutate(across(c(Superior, Michigan, Huron, Erie, Ontario),
                ~ replace_na(.x, "0"))) %>% # convert NA values to 0 for all columns listed
  mutate(across(c(Superior, Michigan, Huron, Erie, Ontario),
                ~ ifelse(.x == "X", 1, 0))) %>%  # convert X to numeric 1 to indicate presence, and any other value, i.e., "extirpated" to 0 for absence
  filter(!if_all(c(Superior, Michigan, Huron, Erie, Ontario),
         ~ .x == 0)) %>% # removes any rows where the species is not present in any of the lakes
  clean_names() # converts all column names to snake_case

# Inspect cleaned data frame
head(df_fish2)
print(df_fish2$species) # noticed few species names were not correctly formatted

# Adjust species names in data frame manually
df_fish2$species <- df_fish2$species %>%
  recode(
    "Chrosomuseos" = "Chrosomus eos",
    "Chrosomuserythrogaster" = "Chrosomus erythrogaster",
    "Chrosomusneogaeus" = "Chrosomus neogaeus", 
    "Sandercanadensis" = "Sander canadensis", 
    "Sandervitreus" = "Sander vitreus", 
    "Lethenteronappendix" = "Lethenteron appendix"
  )

# Inspect the new data frame after names have been adjusted
head(df_fish2)
print(df_fish2$species)

#### _ Exploring and Cleaning Chordate Data from BOLD --------
# Use BOLD API to get nucleotide sequences for all fish - Phylum Chordata
# df_chordata <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Chordata&format=tsv")
# write_tsv(df_chordata, "../data/df_chordata.tsv")

df_chordata <- read_tsv("../data/df_chordata.tsv")

# Inspect chordate data frame
class(df_chordata)
dim(df_chordata)
head(df_chordata)
tail(df_chordata)
str(df_chordata)
names(df_chordata)
table(df_chordata$markercode) # look at what marker codes are available and which ones are the best represented

# Create a subset of the chordata data that includes the relevant information for this analysis
df_chordata2 <- df_chordata %>% 
  dplyr::select(species_name, class_name, markercode, nucleotides) %>% 
  filter(!is.na(species_name)) %>% 
  filter(markercode == "COI-5P") %>% 
  rename(species_name = "species", class_name = "class") %>% # this needs to match the df_fish2 data frame for the left join
  distinct() %>% 
  clean_names()

# Inspect chordate data frame to see how it has changed after filtering
head(df_chordata2)
tail(df_chordata2)
names(df_chordata2)

table(df_chordata2$species)
table(df_chordata2$markercode) # all records should be COI-5P

#### _ Joining Data Frames --------
# Join df_chordata2 with df_fish2
df_fish3 <- inner_join(df_fish2, df_chordata2, by = "species") %>% 
  select(species, class, superior, michigan, huron, erie, ontario, nucleotides)

# Check how many classes are represented and by what proportion
table(df_fish3$class)

# Evaluate if any rows/species were lost in the left join
length(unique(df_fish3$species)) == length(unique(df_fish2$species)) # this evaluates to false, so species were lost
length(unique(df_fish3$species))
length(unique(df_fish2$species)) # there are three more species here, which means that 3 species did not have sequence data

table(df_fish3$species) # since the minimum number of records for a unique species is 1, only keep one record per species to ensure no bias when creating phylogeny 

# Get rid of unnecessary data objects
rm(df_chordata, df_chordata2, df_fish, df_fish2)

#### _ Filtering Nucleotide Sequences --------
# Filtering Nucleotide Sequences to pick Representative Sequence for each Species

##### _ Filter #1: Sequence Length -------- 
# Filter out nucleotide sequences less than 500bp and more than 800bp to ensure we capture COI-5P and not another gene / locus
# count sequence lengths
df_fish_filtered <- df_fish3 %>% 
  group_by(species) %>% 
  mutate(seq_length = nchar(nucleotides)) %>% 
  ungroup()

# visualize the sequence lengths distribution
(hist_sequence_lengths1 <- hist(df_fish_filtered$seq_length))

# filter sequences lss than 500bp and more than 800bp
df_fish_filtered <- df_fish_filtered %>% 
  filter(seq_length > 500 & seq_length < 800)

# visualize the sequence lengths distribution
(hist_sequence_lengths2 <- hist(df_fish_filtered$seq_length))

# remove df_fish3 as it is not needed
rm(df_fish3)

##### _ Filter #2: Gaps and Ambiguous Nucleotides --------
# First add columns to count the number of ambiguous nucleotides and gaps
df_fish_filtered <- df_fish_filtered %>%
  mutate(n_count = str_count(nucleotides, "N"), 
    gap_count = str_count(nucleotides, "-"), 
    nucleotides = as.character(nucleotides)) # convert to character rather than list of characters

# Pick representative: fewest Ns, then longest
df_fish_filtered <- df_fish_filtered %>%
  group_by(species) %>%
  arrange(n_count, gap_count, desc(seq_length)) %>%
  dplyr::slice(1) %>% 
  ungroup()

# Visualize the sequence lengths distribution
hist_sequence_lengths3 <- hist(df_fish_filtered$seq_length)
summary(df_fish_filtered$seq_length)

# Create a 3-panel plot to show the sequence length distribution changes after filtering
png("../fig/01_fig_hist_sequence_lengths_filter.png", width = 2500, height = 700, res = 300)
par(mfrow = c(1, 3),
    mar = c(5, 5, 4, 4),
    cex.lab = 1.0,
    cex.main = 1.0)

# Panel 1 — before filtering
plot(hist_sequence_lengths1,
     main = "A: Before Filtering",
     xlab = "Sequence Length (bp)",
     col = "#3f5a51")

# Panel 2 — after 500–800bp sequence length filter
plot(hist_sequence_lengths2,
     main = "B: After Sequence Length Filter",
     xlab = "Sequence Length (bp)",
     col = "#7FC0BA")

# Panel 3 — after ambiguous base + gap filtering
plot(hist_sequence_lengths3,
     main = "C: After Ambiguous Base and Gap Filter",
     xlab = "Sequence Length (bp)",
     col = "#5aa6a0")

dev.off()

# remove histogram objects 
rm(hist_sequence_lengths1, hist_sequence_lengths2, hist_sequence_lengths3)

# Evaluate the number of entries for each species - should be one for each and ensure there are no missing values in the nucleotide sequences
table(df_fish_filtered$species)
sum(is.na(df_fish_filtered$nucleotides))

#### _ Prepare Data for Analysis -------- 

### Create fasta file with nucleotide sequences and species names
# Prepare intermediate data frame to convert to fasta with header and sequences
df_sequences_fasta <- df_fish_filtered %>%
  transmute(seq.name = gsub(" ", "_", species), seq.text = nucleotides) # make headers for fasta file with underscores to make sure names look proper in MEGA once exported

dat2fasta(dat = df_sequences_fasta, outfile = "../data/phylo/fish_sequences.fasta")

table(df_sequences_fasta$seq.name)

# remove df_sequences_fasta as file is saved and we no longer need it
rm(df_sequences_fasta)

## _ Analyzing Community Phylogenetics -------- 
# NEXT STEPS WERE DONE IN MEGA and IQTREE 2
# 1. align sequences using MUSCLE alignment in MEGA - -400.00 gap open penalty: produced "../data/phylo/MEGA/fish_sequences_aligned_MUSCLE.fas" file
# 2. build ML tree in IQTREE2 - 1000 replicates bootstrap: produced "../data/phylo/IQ-TREE2.fish_sequences_aligned_MUSCLE.fas.contree" file (used below)

### _ Import and Visualize Phylogenetic Tree -------- 
tree_fish <- read.tree("../data/phylo/IQ-TREE2/fish_sequences_aligned_MUSCLE.fas.contree") # this is the tree file created in IQ-TREE 2

# Investigate tree data
str(tree_fish)
summary(tree_fish)

# Adjust tip labels to remove underscores in names before plotting
tree_fish$tip.label <- gsub("_", " ", tree_fish$tip.label)

# Plot phylogenetic tree
png("../fig/02_fig_tree_fish.png", width = 14, height = 14, units = "in", res = 300)
ggtree(tree_fish, layout = "fan") +
  geom_tiplab(size = 2)
dev.off()

### _ Community Structure Analysis --------
# Create community matrix for picante() 
comm_matrix_fish <- df_fish_filtered %>% 
  select(species, superior, michigan, huron, erie, ontario) %>% 
  column_to_rownames("species")

comm_matrix_fish <- as.data.frame(t(comm_matrix_fish))
class(tree_fish)

tree_fish$tip.label

# Check to see if column names in my matrix match the tip labels for my tree
setdiff(colnames(comm_matrix_fish), tree_fish$tip.label)

# Visualize which fish exist in which lakes and make community matrix a data frame to plot 
df_comm_matrix_fish_long <- as_tibble(comm_matrix_fish, rownames = "lake") %>%
  pivot_longer(-lake, names_to = "species", values_to = "presence") %>%
  filter(presence == 1)

fish_species_presence <- ggplot(df_comm_matrix_fish_long, aes(x = species, y = lake)) +
  geom_tile(fill = "#3f5a51") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Great Lake") +
  xlab("Species") +
  ggtitle("Species Presence Across Great Lakes")

ggsave(filename = "03_fig_fish_species_presence.png", plot = fish_species_presence, path = "../fig/", 
       width = 20, height = 8, dpi = 300)

# Remove uneeded objects
rm(df_comm_matrix_fish_long, fish_species_presence)

# Create a phylogenetic distance matrix
fish_dist <- cophenetic.phylo(tree_fish)
head(fish_dist)

# Visualize the distance matrix
png("../fig/04_fig_heatmap_fish_phylo_dist.png", width = 14, height = 14, res = 300, unit = "in")
heatmap(as.matrix(fish_dist), symm = TRUE, col = viridis::viridis(100), margins = c(12, 12))
dev.off()

#### _ Nearest Taxon Index (NTI) Calculation and Visualization --------
# Using Nearest taxon index (NTI) as a measure for phylogenetic structure of communities
# Estimate standardized effect Size of the MNTD on the cophenetic matrix
fish_ses_mntd <- ses.mntd(comm_matrix_fish, fish_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 100)
head(fish_ses_mntd)

capture.output(fish_ses_mntd, file = "../output/01_output_fish_ses_mntd.txt")

# Calculate NTI (z-score from obs-mean/std)
fish_nti <- as.matrix(-1 * ((fish_ses_mntd[,2] - fish_ses_mntd[,3]) / fish_ses_mntd[,4]))
rownames(fish_nti) <- row.names(fish_ses_mntd)
colnames(fish_nti) <- "NTI"

capture.output(fish_nti, file = "../output/02_output_fish_nti.txt")

head(fish_nti)

## Visualize overdispersion / filtering
fish_nti <- as_tibble(fish_nti, rownames = "lake") # convert matrix to tibble data frame object for visualizing 
class(fish_nti)

fish_nti$pattern <- ifelse(fish_nti$NTI > 0, "clustering", "overdispersion") # add pattern of how fish are assembling

lake_cluster_disperse <- ggplot(fish_nti, aes(x = lake, y = NTI, fill = pattern)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.96, linetype = "dashed", colour = "red") + 
  geom_hline(yintercept = -1.96, linetype = "dashed", colour = "red") + 
  annotate("text", x = 0.5, y = 1.96, label = "+1.96", color = "red", hjust = 0, vjust = -0.5) +
  annotate("text", x = 0.5, y = -1.96, label = "-1.96", color = "red", hjust = 0, vjust = 1.5) +
  scale_fill_manual(values = c("clustering" = "#3f5a51", "overdispersion" = "#7FC0BA")) +
  ggtitle("Nearest Taxon Index for Fish Species Across Great Lakes") + 
  ylab("Nearest Taxon Index (z-score)") +
  xlab("Great Lake") +
  labs(fill = NULL) +
  theme_minimal()

ggsave(filename = "05_fig_lake_cluster_disperse.png", plot = lake_cluster_disperse, path = "../fig/")

rm(lake_cluster_disperse)
