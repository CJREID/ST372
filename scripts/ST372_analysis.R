####-----PACKAGES-----####
# Load required packages
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)
library(readr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(pheatmap)
library(ggtree)
library(RColorBrewer)
library(abricateR)
library(plasmidmapR)
library(paletteer)

####-----SETUP-----####
# SET WORKING DIRECTORY
## YOU NEED TO CHANGE THIS TO THE DIRECTORY YOU HAVE CLONED THE GITHUB REPOSITORY TO
wrkdir <- "/Volumes/126451/WORK/Colleagues/Paarthiphan/ST372/data/ST372_ANALYSIS/"
setwd(wrkdir)

# Define 'not in' function for subsetting
'%notin%' <- Negate('%in%')

# Make output directory
if (dir.exists("outputs")){
} else {
  dir.create("outputs")
  dir.create("outputs/data")
  dir.create("outputs/figures")
}

####-----IMPORT DATA-----####
# Get input filenames
infiles <-c(list.files("data/summaries/", pattern = "\\.txt"))

# Read in and name data.frames with whatever comes before .txt
for (file in infiles){
  indir <- c("data/summaries")
  f <- read_delim(paste(indir, file, sep = "/"), delim = "\t", col_names = TRUE, trim_ws = TRUE)
  assign(paste(substr(file, 1, nchar(file)-4), sep = ""), f)
}

# Read in metadata
meta <- read_delim("data/meta/ST372_Metadata_edited.csv", delim = ",") %>% 
  mutate(Source = gsub("Companion Animal", "Canine", Source))

####-----PROCESS ABRICATE DATA-----####
# Load paths to files needed for abricateR
abricate_path <- "data/summaries/genotype.txt"
pointfinder_path <- "data/summaries/pointfinder.txt"
pMLST_data <- "data/summaries/pMLST.txt"

# Provide output names
output_name <- "ST372"
output_dir <- "processed_R"

# Run abricateR to generate gene hits at 90% Id over 95% gene length
abricateR::abricateR(
  abricate_in = abricate_path,
  output = output_name,
  identity = 90,
  length = 95,
  output_directory = output_dir,
  writecsv = FALSE,
  pointfinder_data = pointfinder_path,
  pMLST_data = pMLST_data
)

####-----SEROTYPE DATA-----####
# Read in serotype data
raw_sero <- read_delim("data/summaries/ST372_serotype.tsv")

# Select O_type genes
O_type <- raw_sero %>% 
  select(Name = last_col(), everything(), - Contig, -`Position in contig`, - `Accession number`) %>%
  filter(Database == "O_type") 

# Split to wzx alleles
wzx <- O_type %>% filter(Gene == "wzx")

# Filter sequences that have a single wzx hit
wzx_single <- wzx %>% filter(Name %notin% wzx$Name[duplicated(wzx$Name)])

# Filter sequences that have more than one wzx hit
wzx_dup <- wzx %>% filter(Name %in% wzx$Name[duplicated(wzx$Name)])

# Since all the hits are identical for the different gene names, just select the top one. 
# Most of these are going to be O2/O50 based on their wzy gene too
wzx_dup <- wzx_dup %>% group_by(Name) %>% slice(1)

# Join single and resolved duplicates
wzx <- rbind(wzx_single, wzx_dup)

# Select wzy alleles
wzy <- O_type %>% filter(Gene == "wzy")

# Filter sequences with a single wzy hit
wzy_single <- wzy %>% filter(Name %notin% wzy$Name[duplicated(wzy$Name)])

# Filter sequences with more than one wzy hit (use distinct to remove identical O117 hits for two strains)
wzy_dup <- wzy %>% filter(Name %in% wzy$Name[duplicated(wzy$Name)]) %>% distinct()

# The remaining clashes are between O18 vs O18ac and O2 vs O2/O50 
# I will just go with O18 and O2/O50 as mentioned above for wzx
# Can easily filter to select these
wzy_dup <- wzy_dup %>% filter(Serotype == "O18" | Serotype == "O2/O50" | Serotype == "O117")

# Join single and resolved duplicates
wzy <- rbind(wzy_single, wzy_dup)

# Rejoin wzy and wzx alleles and sort by Name
O_type <- rbind(wzx, wzy) %>% arrange(Name)

## Now we have samples with a wzy and a wzx hit or one of the two. 
## Obviously the O-type for the single hit ones is resolved but the wzx/wzy hits 
## need to be compared for those samples with both alleles 

# Get names of smaples with two hits
two_alleles <- O_type %>% count(Name) %>% filter(n >1) %>% pull(Name)

# Filter O_type double hits that are in agreement and define their O_type
O_type_a <- O_type %>% 
  filter(Name %in% two_alleles) %>%
  group_by(Name, Serotype) %>% summarise(sero_count = n()) %>%
  filter(sero_count == 2) %>%
  select(Name, O_type = Serotype)

# Filter O_type to double hit samples with conflicting types
conflicting_O <- O_type %>%
  filter(Name %in% two_alleles) %>% 
  group_by(Name, Serotype) %>% 
  summarise(sero_count = n()) %>% 
  filter(sero_count == 1)

# All these are O18 vs O18ac so I will select O18 for simplicity, because who the hell knows what the difference is?
O_type_b <- conflicting_O %>% select(Name, O_type = Serotype) %>% filter(O_type == "O18")

# Filter single hit samples
O_type_c <- O_type %>% filter(Name %notin% two_alleles) %>% select(Name, O_type = Serotype)

# Join them all together
O_type <- rbind(O_type_a, O_type_b, O_type_c)

# Get H types
H_type <- raw_sero %>% select(Name = last_col(), Database, Serotype)%>%
  filter(Database == "H_type") %>% select(Name, H_type=Serotype)

# Join O and H to complete serotype data
serotype <- left_join(H_type, O_type) %>% 
  mutate(O_type = case_when(is.na(O_type) ~ "ONT", TRUE ~ O_type)) %>%
  select(Name, O_type, H_type) %>%
  mutate(Serotype = paste(O_type, H_type, sep = ":"))


# Split into types with 10 or more, and less than 10 representatives
major_sero <- serotype %>% select(Name, Serotype) %>% group_by(Serotype) %>% filter(n() >= 10)
minor_sero <- serotype %>% select(Name, Serotype) %>% group_by(Serotype) %>% filter(n() < 10)

# Designate types with less than 10 as "Other"
minor_sero$Serotype <- "Other OH"

# Stick them back together
serotype_simple <- rbind(major_sero, minor_sero) %>% select(Name, Serotype_simple = Serotype)

####-----PROCESS FIMH DATA-----####
# Read in raw fimH data and select the useful columns
fimH_raw <- read_csv("data/summaries/ST372.fimH.csv") %>% 
  mutate(Name = gsub(".fasta", "", `#FILE`)) %>%
  select(Name, fimH = GENE, Coverage = `%COVERAGE`, ID = `%IDENTITY`)

# There are more rows than samples so there must be sequences with duplicate hits
# Get subset that only has single hits
fimH_single <- fimH_raw %>% filter(Name %notin% fimH_raw$Name[duplicated(fimH_raw$Name)]) %>% distinct()

# Get the subset with more than one hit --> they're all identical alleles hitting
# different ranges on the query so can be made single hits with distinct()
fimH_dup <- fimH_raw %>% filter(Name %in% fimH_raw$Name[duplicated(fimH_raw$Name)]) %>% distinct()

# Join them back together 
fimH <- rbind(fimH_single, fimH_dup) %>% arrange(Name) %>% select(Name, fimH)

# Split into types with 10 or more, and less than 10 representatives
major_fimH <- fimH %>% select(Name, fimH) %>% group_by(fimH) %>% filter(n() >= 10)
minor_fimH <- fimH %>% select(Name, fimH) %>% group_by(fimH) %>% filter(n() < 10)

# Designate types with less than 10 as "Other"
minor_fimH$fimH <- "Other fimH"

# Re-join them
fimH_simple <- rbind(major_fimH, minor_fimH) %>% arrange(Name) %>% select(Name, fimH_simple = fimH)

####-----PROCESS pUTI89 ALIGNMENT-----####
# Get tree path and abricate alignment data
path_to_tree <- "data/iqtree/core_gene_alignment.OGremoved.contree"
path_to_abricate <- "data/summaries/ST372.pUTI89.tab"
plasrefname <- "pUTI89"

# Define reference length so we can calculate %ID later
pUTI89_ref_length <- read_delim(path_to_abricate, n_max = 1) %>% 
  mutate(COVERAGE = as.integer(gsub("^.*\\/", "", COVERAGE))) %>% pull(COVERAGE)

# Run plasmidmapR to generate binned hit table for each ST372 sequence against pUTI89 sequence
plasmid_mapR(path_to_abricate = path_to_abricate,
             plasmid_reference_name = plasrefname,
             output_directory = NULL,
             min_hit_id = 90,
             min_hit_length = 0.5,
             writecsv = FALSE,
             path_to_tree = path_to_tree)

# Generate dataframe with a summary statistic for how many bins were covered at >90%ID
pUTI89_ID <- as.data.frame(rowSums(pUTI89_binned_hits))

# Add names column
pUTI89_ID$Name <- rownames(pUTI89_binned_hits)

# Rename columns
colnames(pUTI89_ID) <- c("pUTI89_ID","Name")

# Convert value to a percentage
pUTI89_ID$pUTI89_ID <- round((pUTI89_ID$pUTI89_ID/pUTI89_ref_length) * 100)

####-----TREE-----####
# Read in the tree file
tree <- read.tree(file = path_to_tree)

####-----BAPS-----####
# Read in BAP cluster data
bap_gene <- read_csv("data/fastbaps/ST372.fastbaps.csv") %>% select(Name = Isolates, everything())

# Extract names
strain_names <- bap_gene$Name

# Re-assign names to rownames
rownames(bap_gene) <- strain_names

# Generate data frame for plotting panels for BAP clustering of each alignment
plot.df.gene <- data.frame(id = rownames(bap_gene), 
                           fastbaps1 = bap_gene$`Level 1`, 
                           stringsAsFactors = FALSE)

# Create a character column of the values to make plotting easier
plot.df.gene$gene_bap <- as.character(plot.df.gene$fastbaps1)

# Add BAP prefix to the cluster numbers
plot.df.gene$gene_bap <- paste0("BAP", plot.df.gene$gene_bap)

# Create ggtrtee object from ST127 phylogenetic tree
bap.tree <- ggtree(tree, branch.length = "none") 

# Create facet plots showing fastbaps grouping for each phylogeny
facet_gene <- facet_plot(bap.tree, 
                         panel = "fastbaps",
                         data = plot.df.gene,
                         geom = geom_tile, 
                         aes(x = fastbaps1, 
                             fill = gene_bap))

# Reorder plot.df.gene so that we can rename the clusters with sequential letters based on the tree
d <- fortify(tree)
dd <- subset(d, isTip)
colorder <- dd$label[order(dd$y, decreasing=TRUE)]
colordnum1 <- match(colorder, plot.df.gene$id)
bap_rename <- plot.df.gene[colordnum1,]

# Rename the clusters to letters that follow the tree structure
bap_rename <- bap_rename %>% 
  mutate(Cluster = case_when(gene_bap == "BAP1" ~ "A",
                             gene_bap == "BAP13" ~ "B",
                             gene_bap == "BAP11" ~ "C",
                             gene_bap == "BAP14" ~ "D",
                             gene_bap == "BAP10" ~ "E",
                             gene_bap == "BAP3" ~ "F",
                             gene_bap == "BAP9" ~ "G",
                             gene_bap == "BAP7" ~ "H",
                             gene_bap == "BAP2" ~ "I",
                             gene_bap == "BAP6" ~ "J",
                             gene_bap == "BAP4" ~ "K",
                             gene_bap == "BAP5" ~ "L",
                             gene_bap == "BAP8" ~ "M"))

# Process data into casted format for Scoary analysis prior to excluding the outgroup strain
baps.scoary <- bap_rename

# Create value column for casting
baps.scoary$value <- rep(1, nrow(bap_rename))

# Select relevant columns
baps.scoary <- baps.scoary %>% select(Name = id, Cluster, value)

# Cast into wide format
baps.scoary <- baps.scoary %>% reshape2::dcast(Name ~ Cluster)

# Recode NAs to zeroes
baps.scoary[is.na(baps.scoary)] <- 0

# Write file for Scoary
# write_csv(baps.scoary, "outputs/data/ST372.BAPS.scoary.csv")

# Add BAP cluster data to metadata
meta <- left_join(meta, bap_rename, by = c("Name" = "id")) %>% select(everything(), -fastbaps1, Cluster)

####-----PROCESS METADATA-----####
# Select ColV data to add to meta
colv <- ST372_simple_summary_N90L95 %>% select(Name = name, ColV, IncF_RST)

# Add ST data - all ST372 but add to be sure, to be sure
mlst <- read_delim("data/summaries/mlst.txt")

# Confirm all ST372
meta <- left_join(meta, mlst, by = c("Name" = "name")) %>% mutate(ST = str_extract(ST, "^[0-9]{3}")) %>% select(-scheme)

# Join plasmid meta and drop STs
meta <- left_join(meta, colv)

# Join serotype data
meta <- left_join(meta, serotype)

# Join simplified serotype data
meta <- left_join(meta, serotype_simple)

# Join fimH data
meta <- left_join(meta,fimH)

# Join simplified fimH data
meta <- left_join(meta,fimH_simple)

# Categorise ColV
meta <- meta %>% mutate(ColV = case_when(ColV == "0" ~ "No", ColV == "1" ~ "Yes"))

# Categorise seqeunces that are F29:A-:B10 carriers or >90% coverage of pUTI89 as pUTI89+
meta <- left_join(meta, pUTI89_ID) %>%
  mutate(pUTI89 = case_when(pUTI89_ID >= 90 | IncF_RST == "F29:A-:B10" ~ "Yes", TRUE ~ "No")) %>%
  select(-pUTI89_ID)

# Extract working names for collection
working_names <- as.vector(meta$Name)

# Create 'geno_meta' df with all meta and gene screening data
geno_meta <- left_join(meta %>% select(Name, Source, Year, Continent, Country, ColV, pUTI89, Cluster), 
                       ST372_simple_summary_N90L95 %>% select(-ColV), by = c("Name" = "name"))

# Define gene columns as those that are integers
gene_cols <- names(geno_meta %>% select(where(is.integer)))

# Recode multiple hits as a single hit
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(2, 1, .x)))
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(3, 1, .x)))
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(4, 1, .x)))

# Convert back to integer (gsub makes everything character class)
geno_meta <- geno_meta %>% mutate(across(all_of(gene_cols), as.integer))

# Extract metadata column names for later use with the geno_meta dataframe
meta_cols <- c("Name", "Source", "Year", "Continent", "Country", "ColV", "pUTI89", "Cluster")

####-----TREE METADATA-----####
# Extract just what we need for tree visualisation
tree_meta <- meta %>% 
  select(Name,
         Cluster,
         Source,
         Year,
         Continent, 
         Country,
         ColV,
         IncF_RST,
         pUTI89,
         fimH = fimH_simple,
         Serotype = Serotype_simple)

# Get names
tree_names <- tree_meta$Name

# Refine tree metadata
tree_meta <- tree_meta %>% select(Cluster, Source, Continent, Serotype, fimH)

# Assign names as row names (for gheatmap function)
rownames(tree_meta) <- tree_names

####-----DEFINE COLOURS-----####
# pUTI89 colours
puti_vars <- unique(meta$pUTI89)
puti_clrs <- c("#38c9b1", "white")
names(puti_clrs) <- rev(puti_vars)

# Colours for BAP clusters
bap_vars <- unique(meta$Cluster)
bap_clrs <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(meta$Cluster)))
names(bap_clrs) <- bap_vars

# Colours for isolate source
source_vars <- unique(meta$Source)
source_clrs <- colorRampPalette(brewer.pal(length(unique(meta$Source)), "Set1"))(length(unique(meta$Source)))
names(source_clrs) <- sort(source_vars)
source_clrs["Canine"] <- "#377EB8CC"
source_clrs["Human"] <- "#984EA3E6"

# Colours for continents
continent_vars <- unique(meta$Continent)
continent_clrs <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(meta$Continent)))
set.seed(8)
continent_clrs <- sample(continent_clrs)
names(continent_clrs) <- sort(continent_vars)

# Colours for year
year_vars <- unique(meta$Year)
year_clrs <- colorRampPalette(paletteer_dynamic("cartography::blue.pal", 20))(length(unique(meta$Year)))
names(year_clrs) <- sort(year_vars)
year_clrs["ND"] <- "#dbdbdb"

# Colours for serotypes
OH_vars <- unique(meta$Serotype_simple)
OH_clrs <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(meta$Serotype_simple)))
names(OH_clrs) <- OH_vars

# Colours for fimH
fim_vars <- unique(meta$fimH_simple)
fim_clrs <- colorRampPalette(brewer.pal(5, "Set1"))(length(unique(meta$fimH_simple)))
names(fim_clrs) <- sort(fim_vars)

# Combine colours for gheatmap
tree_vars <- c(source_clrs, continent_clrs, bap_clrs, OH_clrs, fim_clrs)

####-----FIGURE 1 METADATA-----####
# Plot isolate continents of origin stratified by source
fig1a <- ggplot(meta, aes(Continent)) +
  geom_bar(aes(fill = Source))+
  geom_text(stat='count', aes(label=..count..), vjust=-0.3, size = 3) +
  scale_fill_manual(name = "Source", values = source_clrs) +
  scale_x_discrete(name = "Continent", labels = function(x) str_wrap(x, width = 10))+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 270), n.breaks = 8) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 9, vjust = , face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6))

# Plot clusters stratified by source
fig1b <- ggplot(meta, aes(Cluster)) +
  geom_bar(aes(fill = Source))+
  geom_text(stat='count', aes(label=..count..), vjust=-.3, size = 3) +
  scale_fill_manual(name = "Source", values = source_clrs) +
  scale_x_discrete(name = "Cluster", labels = function(x) str_wrap(x, width = 10))+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 220), n.breaks = 8) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6))

# Plot sources stratified by serotypes
fig1c <- ggplot(meta, aes(Source)) +
  geom_bar(aes(fill = Serotype_simple))+
  geom_text(stat='count', aes(label=..count..), vjust=-.3, size = 3) +
  scale_fill_manual(name = "Serotype", values = OH_clrs) +
  scale_x_discrete(name = "Source", labels = function(x) str_wrap(x, width = 10))+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 330), n.breaks = 8) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 9, vjust = , face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6)) 

# Plot clusters stratified by serotypes
fig1d <- ggplot(meta, aes(Cluster)) +
  geom_bar(aes(fill = Serotype_simple))+
  geom_text(stat='count', aes(label=..count..), vjust=-.3, size = 3) +
  scale_fill_manual(name = "Serotype", values = OH_clrs) +
  scale_x_discrete(name = "Cluster", labels = function(x) str_wrap(x, width = 10))+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 220), n.breaks = 8) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6))

# Plot source stratified by fimH alleles
fig1e <- ggplot(meta, aes(Source)) +
  geom_bar(aes(fill = fimH_simple))+
  geom_text(stat='count', aes(label=..count..), vjust=-.3, size = 3) +
  scale_fill_manual(name = "fimH", values = fim_clrs) +
  scale_x_discrete(name = "Source", labels = function(x) str_wrap(x, width = 10))+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 330), n.breaks = 8) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 9, vjust = , face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6)) 

# Plot clusters stratified by fimH alleles
fig1f <- ggplot(meta, aes(Cluster)) +
  geom_bar(aes(fill = fimH_simple))+
  geom_text(stat='count', aes(label=..count..), vjust=-.3, size = 3) +
  scale_fill_manual(name = "fimH", values = fim_clrs) +
  scale_x_discrete(name = "Cluster", labels = function(x) str_wrap(x, width = 10))+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 220), n.breaks = 8) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6))

# Combine a and b for multiplot
fig1ab <- ggarrange(plotlist= list(fig1a, fig1b),
                    ncol =2, labels =c("a)", "b)"), 
                    font.label = list(size =10),
                    common.legend = TRUE,
                    legend = "right",
                    align = "hv")

# Combine c and d for multiplot
fig1cd <- ggarrange(plotlist= list(fig1c, fig1d),
                    ncol =2, labels =c("c)", "d)"),
                    font.label = list(size =10),
                    legend = "right", common.legend = TRUE, align = "hv")

# Combine e and f for multiplot
fig1ef <- ggarrange(plotlist= list(fig1e, fig1f),
                    ncol =2, labels =c("e)", "f)"),
                    font.label = list(size =10),
                    legend = "right", common.legend = TRUE, align = "hv")

# Plot all together for figure 1
fig1 <- ggarrange(plotlist = list(fig1ab, fig1cd, fig1ef),
                  nrow =3, align = "hv")

# Save as Figure 1
ggsave("Figure1.metadata.pdf",
       fig1, 
       path = "outputs/figures/", 
       device = "pdf", 
       width= 297, 
       height = 180, 
       unit ="mm", 
       dpi = 300)

####-----FIGURE 2 TREE-----####
# Create ggtree object from phylogenetic tree
tree.heat <- rotate_tree(ggtree(tree, layout = "fan", open.angle = 7, size = .2),91)

# Plot tree and metadata
fig2 <- gheatmap(tree.heat, 
                 tree_meta, 
                 width = .2,
                 font.size = 1.7,
                 colnames_offset_x = ,
                 colnames_offset_y = 3.6,
                 colnames_position = "top",
                 colnames_angle = ,
                 hjust = 0.5,
                 color = NULL) +
  scale_fill_manual(name = "Data", values = tree_vars) + 
  theme(legend.position = "none",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10))

# Save as Figure 2
ggsave("Figure2.tree.pdf",
       fig2, 
       path = "outputs/figures/", 
       device = "pdf", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 600)

####-----FIGURE 3 SNP ANALYSIS-----####
# Read in pairwise SNP-dists file
core_snps <- read_csv(file = "data/snps/CGA.snp-dists.csv") %>% rename(Name = ...1)

# Extract names
rows <- core_snps$Name

# Melt into long format
melt.core <- reshape2::melt(core_snps)

# Join metadata to first isolate in each pair
melt.core.a <- left_join(melt.core, meta, by=c("Name")) %>% 
  select(Name, variable, value, Source, Continent, Cluster) %>% 
  rename(Name1 = Name,
         Name = variable,
         Source1 = Source,
         Continent1 = Continent,
         Clustera = Cluster)

# Join metadata to second isolate in each pair, then filter to equal to or less than 30 SNPs. Exclude self comparisons by filtering out 0 SNPS
# Filter out identical source pairs
# Filter Source 1 to Human and Companion animal and Source two to anything but Companion animals
# This results in 77 unique pairs of human or companion animal vs all comparisons
melt.core.unique <- left_join(melt.core.a, meta, by=c("Name")) %>%
  select(Name1, Name, Source1, Source, Continent1, Continent,  Clustera, Cluster, value)%>% 
  rename(Name2 = Name, Source2 = Source, Continent2 = Continent,Clusterb = Cluster, `Core SNPs` = value) %>% 
  filter(`Source1` != `Source2`, `Core SNPs` <=30, `Core SNPs` != 0, Source1 == "Human" | Source1 == "Canine", Source2 != "Canine")

# Extract the countries for visualisation later
# country1 <-  melt.core.unique %>% filter(`Core SNPs` <= 30) %>% pull(Country1) %>% unique()
# country2 <-  melt.core.unique %>% filter(`Core SNPs` <= 30) %>% pull(Country2) %>% unique()
# snp_countries <- unique(c(country1, country2))

# Extract annotation row data for pheatmap
anno_rows.unique <- melt.core.unique %>% 
  filter(`Core SNPs` <= 30) %>% 
  select(Name1, 
         Source = Source1,
         Continent = Continent1,
         Cluster = Clustera) %>%
  distinct()

# Convert Name column to rownames
rownem.unique <- anno_rows.unique$Name1
anno_rows.unique <- anno_rows.unique %>% select(-Name1)
rownames(anno_rows.unique) <- rownem.unique

# As above except for columns
anno_cols_snp.unique <- melt.core.unique %>% 
  filter(`Core SNPs` <= 30) %>%
  select(Name2, 
         Source = Source2,
         Continent = Continent2,
         Cluster = Clusterb) %>% 
  distinct()

# Convert name to rownames
colnem.unique <-anno_cols_snp.unique$Name2
anno_cols_snp.unique <- anno_cols_snp.unique %>% select(-Name2)
rownames(anno_cols_snp.unique) <- colnem.unique

# Create DF of paired names and SNP counts
pairs.unique <- melt.core.unique %>% select(Name1, `Core SNPs`, Name2) %>% filter(`Core SNPs` <= 30)

# Cast to SNP matrix
paircast.unique <- reshape2::dcast(pairs.unique, Name1 ~ Name2, value.var = "Core SNPs")

# Move Name column to rownames and convert to matrices
SNP_rows.unique <- paircast.unique$Name1
paircast_mat.unique <- as.matrix(paircast.unique %>% select(-Name1))
rownames(paircast_mat.unique) <- SNP_rows.unique

# Set NAs outside the SNP range so they can be ignored
paircast_mat.unique[is.na(paircast_mat.unique)] <- 100

# Make quick heatmap with clustering to group the most closely related strains
out.unique <- pheatmap(mat = paircast_mat.unique,
                       cluster_cols = TRUE,
                       cluster_rows = TRUE)

# Extract row and column orders from the above heatmap
roword.unique <- rownames(paircast_mat.unique[out.unique$tree_row[["order"]],])
colord.unique <- colnames(paircast_mat.unique[,out.unique$tree_col[["order"]]])

# Reorder matrix with the order extracted above
paircast_mat.unique <- paircast_mat.unique[roword.unique,]
paircast_mat.unique <- paircast_mat.unique[,colord.unique]

# Reset 100 values as NAs
paircast_mat.unique[paircast_mat.unique == 100] <- NA

# Define country colours for heatmap annotation
country_clrs <- colorRampPalette(brewer.pal(8, "Set1"))(length(snp_countries))
names(country_clrs) <- snp_countries

# Combine all colours for heatmap annotation
anno_colors <- list(Source = source_clrs,
                    Continent = continent_clrs, Cluster = bap_clrs)

# Create Figure 3
fig3 <- pheatmap(mat = paircast_mat.unique,
                     cluster_cols = FALSE,
                     cluster_rows = FALSE,
                     annotation_col = anno_cols_snp.unique,
                     annotation_row = anno_rows.unique,
                     color = colorRampPalette(brewer.pal(9, "BuPu"))(14),
                     breaks = c(seq(from = 0, to = 30, by = 2)),
                     fontsize_col = 8,
                     fontsize_row =8,
                     show_rownames = TRUE,
                     show_colnames = TRUE,
                     annotation_colors = anno_colors,
                     annotation_names_col = TRUE,
                     annotation_names_row = TRUE,
                     legend = FALSE,
                     annotation_legend = FALSE,
                     border_color = "black",
                     na.cols = "grey")

# Save figure 3
ggsave("Figure3.SNP.heatmap.pdf", 
       fig3, 
       path = "outputs/figures/", 
       device = "pdf", 
       width = 297,
       height = 297,
       unit ="mm", 
       dpi = 300)

####-----SCOARY CLUSTER PROCESSING-----####
# Generate file list of Scoary results for BAP groups
filelist <- list.files(path="data/scoary/cluster", pattern="*.csv", full.names = TRUE)

# Read in thefiles that actually have data - used file size as a proxy for this
for (f in 1:length(filelist)){
  if(file.size(filelist[f]) > 229) {
    assign(paste0("cluster_", str_extract(filelist[f], "[A-Z]{1}"), "_scoary"), read_csv(filelist[f]))
  }
}

# Combine into list
bap_scoary_list <- mget(ls(pattern = "cluster_.*_scoary"))

# Processing loop for each file
for (f in 1:length(bap_scoary_list)){
  data <- bap_scoary_list[[f]]
  name <- str_extract(names(bap_scoary_list[f]), "[A-Z]{1}")
  # Filter hypotheticals and combine gene and non-unique gene name to create a unique name that splits paralogs
  data <- data %>% 
    mutate(Cluster = rep(name, nrow(data)), 
           Gene_Unique = paste(`Non-unique Gene name`, Gene, sep = "_")) %>% 
    filter(!grepl("hypothetical", Annotation), Benjamini_H_p < 1E-20) %>% 
    select(Gene_Unique, everything())
  
  # Get rid of NAs
  data$Gene_Unique <- gsub("NA_", "", data$Gene_Unique)
  
  # Split into over and under-represented ColV clade genes based on the Odds Ratio
  over_rep <- data %>% filter(Odds_ratio > 1 | Odds_ratio == "inf") %>%
    # slice(0:20) %>%
    dplyr::select(-Sensitivity,-Specificity,-Bonferroni_p) %>%
    dplyr::rename(
      Pos_present = Number_pos_present_in,
      Neg_present = Number_neg_present_in,
      Pos_absent = Number_pos_not_present_in,
      Neg_absent = Number_neg_not_present_in
    ) 
  
  under_rep <- data %>% filter(Odds_ratio < 1) %>%
    # slice(0:20) %>%
    dplyr::select(-Sensitivity,-Specificity,-Bonferroni_p) %>%
    dplyr::rename(
      Pos_present = Number_pos_present_in,
      Neg_present = Number_neg_present_in,
      Pos_absent = Number_pos_not_present_in,
      Neg_absent = Number_neg_not_present_in
    ) 
  # Assign names to the outputs
  assign(paste(str_extract(names(bap_scoary_list[f]), "cluster_[A-Z]{1}"), "over", sep = "_"), over_rep)
  assign(paste(str_extract(names(bap_scoary_list[f]), "cluster_[A-Z]{1}"), "under", sep = "_"), under_rep)
  
  # Remove unnecessary objects
  rm(data)
  rm(over_rep)
  rm(under_rep)
  rm(name)
}

# Get a list of the data frames from the above loop
bap_output_list <- mget(ls(pattern = "cluster_.*_[over|under]"))

# Create empty list for all significant BAP genes
all_bap_list <- list()

# Process all significant BAP genes into one dataframe
for (f in 1:length(bap_output_list)){
  if(nrow(bap_output_list[[f]]) != 0){
    all_bap_list <- append(all_bap_list, bap_output_list[f])
    all_bap <- bind_rows(all_bap_list)
  }
}

####-----ROARY PROCESSING-----####
# Read in Roary gene_presence_absence
roary_raw_meta <- read_csv("data/roary/ST372.roary.meta.csv", 
                           col_names = TRUE,
                           col_types = cols(
                             .default = col_character(),
                             `Non-unique Gene name` = col_character(),
                             `No. isolates` = col_double(),
                             `No. sequences` = col_double()), 
                           skip_empty_rows = FALSE, guess_max = Inf)

# Read in binary gene presence absence
roary_tab <- read.table("data/roary/gene_presence_absence.Rtab", header =TRUE, check.names = FALSE)

# Join metadata and binary table
roary_filter <- left_join(roary_raw_meta, roary_tab)

# Create percentage carriage column to filter on
roary_filter <- roary_filter %>% 
  mutate(Carriage = `No. isolates`/408*100, # the number 408 needs to be changed based on your collection if reusing this script
         Gene_Unique = paste(`Non-unique Gene name`, Gene, sep = "_")) %>% 
  select(Gene_Unique, 3:4, Carriage, everything()) %>% rename_with(~ gsub(pattern = ".out", "", .x))

# Get rid of NAs
roary_filter$Gene_Unique <- gsub("NA_", "", roary_filter$Gene_Unique)

# Get gene presence absence from full Roary table
roary_filter2 <- roary_filter %>% select(Gene_Unique, where(is.integer))

####-----COMBINE ROARY AND SCOARY DATA-----####
####--CLUSTER M--####
# Select cluster M significant genes from fastBAPS data
clusterM_map_A <- all_bap %>% filter(Cluster == "M") %>% select(Gene_Unique, Annotation, Odds_ratio, Benjamini_H_p) %>% filter(Odds_ratio >=1)

# Pull out the presence/absence of the significant genes
clusterM_map <- left_join(clusterM_map_A, roary_filter2) %>% select(Gene = Gene_Unique, everything(), -Annotation, -Benjamini_H_p, -Odds_ratio)

# Get the gene names for rownames later
clusterMgenes <- c(clusterM_map$Gene)

# Transpose so genes are columns and de-select sample names
clusterM_map <- t(clusterM_map) %>% as.data.frame %>% slice(2:n())

# Replace Vx colnames with gene names
names(clusterM_map) <- clusterMgenes

# Sort names alphabetically
clusterM_map <- 
  clusterM_map %>% select(sort(names(.)))

# Change 1's and zeroes to PresentM/Absent - the 'M' allows us to colour the heatmap 
# differentially later on
clusterM_map[clusterM_map > 0] <- "PresentM"
clusterM_map[clusterM_map ==0] <- "Absent"

# Replace the long 'group_nnnn' names with 'b' to denote alternative allele of a gene
names(clusterM_map) <- gsub("_group.*", "_b", names(clusterM_map))

# Create metadata for downstream heatmap
clusterM_meta <- left_join(meta %>% 
                             select(Name, Cluster, Source, Serotype = Serotype_simple), 
                           clusterM_map %>% 
                             mutate(Name = rownames(clusterM_map)) %>% 
                             select(Name, everything()))

# Create another tree object for plotting
scoary_tree <- ggtree(tree,branch.length = "none", size=.1)

####--CLUSTER G--####
# Select cluster G significant genes from fastBAPS data
clusterG_map_A <- all_bap %>% filter(Cluster == "G") %>% select(Gene_Unique, Annotation, Odds_ratio, Benjamini_H_p) %>% filter(Odds_ratio >=1)

# Pull out the presence/absence of the significant genes
clusterG_map <- left_join(clusterG_map_A, roary_filter2) %>% select(Gene = Gene_Unique, everything(), -Annotation, -Benjamini_H_p, -Odds_ratio)

# Get the gene names for rownames later
clusterGgenes <- c(clusterG_map$Gene)

# Transpose so genes are columns and de-select sample names
clusterG_map <- t(clusterG_map) %>% as.data.frame %>% slice(2:n())

# Replace Vx colnames with gene names
names(clusterG_map) <- clusterGgenes

# Sort colnames
clusterG_map <- 
  clusterG_map %>% select(sort(names(.)))

# Replace 1's and 0's with PresentG/Absent for heatmap
clusterG_map[clusterG_map > 0] <- "PresentG"
clusterG_map[clusterG_map ==0] <- "Absent"

# Replace the long 'group_nnnn' names with 'b' to denote alternative allele of a gene
names(clusterG_map) <- gsub("_group.*", "_b", names(clusterG_map))

# Create metadata for heatmap
clusterG_meta <- left_join(meta %>% 
                             select(Name), 
                           clusterG_map %>% 
                             mutate(Name = rownames(clusterG_map)) %>% 
                             select(Name, everything()))

####--CLUSTER L--####
# Select cluster L significant genes from fastBAPS data
clusterL_map_A <- all_bap %>% filter(Cluster == "L") %>% select(Gene_Unique, Annotation, Odds_ratio, Benjamini_H_p) %>% filter(Odds_ratio >=1)

# Pull out the presence/absence of the significant genes
clusterL_map <- left_join(clusterL_map_A, roary_filter2) %>% select(Gene = Gene_Unique, everything(), -Annotation, -Benjamini_H_p, -Odds_ratio)

# Get the gene names for rownames later
clusterLgenes <- c(clusterL_map$Gene)

# Transpose so genes are columns and de-select sample names
clusterL_map <- t(clusterL_map) %>% as.data.frame %>% slice(2:n())

# Replace Vx colnames with gene names
names(clusterL_map) <- clusterLgenes

# Sort column names
clusterL_map <- 
  clusterL_map %>% select(sort(names(.)))

# Replace 1's and 0's with PresentL/Absent for heatmap
clusterL_map[clusterL_map > 0] <- "PresentL"
clusterL_map[clusterL_map ==0] <- "Absent"

# Replace the long 'group_nnnn' names with 'b' to denote alternative allele of a gene
names(clusterL_map) <- gsub("_group.*", "_b", names(clusterL_map))

# Create metadata
clusterL_meta <- left_join(meta %>% 
                             select(Name), 
                           clusterL_map %>% 
                             mutate(Name = rownames(clusterL_map)) %>% 
                             select(Name, everything()))

####--CLUSTER J--####
# Select cluster J significant genes from fastBAPS data
clusterJ_map_A <- all_bap %>% filter(Cluster == "J") %>% select(Gene_Unique, Annotation, Odds_ratio, Benjamini_H_p) %>% filter(Odds_ratio >=1)

# Pull out the presence/absence of the significant genes
clusterJ_map <- left_join(clusterJ_map_A, roary_filter2) %>% select(Gene = Gene_Unique, everything(), -Annotation, -Benjamini_H_p, -Odds_ratio)

# Get the gene names for rownames later
clusterJgenes <- c(clusterJ_map$Gene)

# Transpose so genes are columns and de-select sample names
clusterJ_map <- t(clusterJ_map) %>% as.data.frame %>% slice(2:n())

# Replace Vx colnames with gene names
names(clusterJ_map) <- clusterJgenes

# Sort colnames
clusterJ_map <- 
  clusterJ_map %>% select(sort(names(.)))

# Replace 1's and 0's with PresentJ/Absent for heatmap
clusterJ_map[clusterJ_map > 0] <- "PresentJ"
clusterJ_map[clusterJ_map ==0] <- "Absent"

# Replace the long 'group_nnnn' names with 'b' to denote alternative allele of a gene
names(clusterJ_map) <- gsub("_group.*", "_b", names(clusterJ_map))

# Create meta for heatmap
clusterJ_meta <- left_join(meta %>% 
                             select(Name), 
                           clusterJ_map %>% 
                             mutate(Name = rownames(clusterJ_map)) %>% 
                             select(Name, everything()))

####-----FIGURE 4 SCOARY GENES WITH PHYLOGENY-------####
# Create colours for presence+cluster association and absence
heat_clrs_MGLJ <- c("PresentM" = "#8DD3C7", 
                    "PresentG" = "#EB8E8B", 
                    "PresentL" = "#A9A0B2", 
                    "PresentJ" = "#D8C965",
                    "Absent" = "#ededed")

# Join the metadata and genes for each cluster
MGLJ_tree_meta <- left_join(clusterM_meta, clusterG_meta)
MGLJ_tree_meta <- left_join(MGLJ_tree_meta, clusterL_meta)
MGLJ_tree_meta <- left_join(MGLJ_tree_meta, clusterJ_meta)

# Extract names
tree_names <- MGLJ_tree_meta$Name

# Remove Name column
MGLJ_tree_meta <- MGLJ_tree_meta %>%select(-Name)

# Make rownames sequence names
rownames(MGLJ_tree_meta) <- tree_names

# Plot gene tree with over-represented genes for clusters M, G, L and J
gheatmap(scoary_tree,
         data = MGLJ_tree_meta,
         font.size = 1.2,
         hjust = 1,
         colnames = TRUE,
         colnames_position = "bottom",
         colnames_angle = 90,
         colnames_offset_y = ,
         colnames_offset_x = ,
         width = 18,
         offset = -3,
         color = NULL) + 
  scale_fill_manual(values = c(source_clrs, bap_clrs, OH_clrs, heat_clrs_MGLJ)) +
  theme(legend.position = "none") + ylim(-20, NA)

# Scale factor for plotting
scalefactor <- 2

# Copy plot
dev.copy(png, res = 300*scalefactor, height=2480, width=3508, units = , 'outputs/figures/Figure4.scoary.heatmap.png')

# Save plot and dev off
dev.off()

####-----ISLANDVIEWER PROCESSING----####
#### pdu_G island ####
# Get tree path and abricate alignment data
path_to_tree <- "data/iqtree/core_gene_alignment.OGremoved.contree"
path_to_abricate <- "data/islandviewer/abricate_screening/ST372.pdu_G.txt"
plasrefname <- "pdu_G"

# Define reference length so we can calculate %ID later
pdu_G_ref_len <- read_delim(path_to_abricate, n_max = 1) %>% 
  mutate(COVERAGE = as.integer(gsub("^.*\\/", "", COVERAGE))) %>% pull(COVERAGE)

# Run plasmidmapR to generate binned hit table for each ST372 sequence against pUTI89 sequence
plasmid_mapR(path_to_abricate = path_to_abricate,
             plasmid_reference_name = plasrefname,
             output_directory = NULL,
             min_hit_id = 90,
             min_hit_length = 0.5,
             writecsv = FALSE,
             path_to_tree = path_to_tree)

# Generate dataframe with a summary statistic for how many bins were covered at >90%ID
pdu_G_ID <- as.data.frame(rowSums(pdu_G_binned_hits))

# Add names column
pdu_G_ID$Name <- rownames(pdu_G_binned_hits)

# Rename columns
colnames(pdu_G_ID) <- c("pdu_G_ID","Name")

# Convert value to a percentage
pdu_G_ID$pdu_G_ID <- round((pdu_G_ID$pdu_G_ID/pdu_G_ref_len) * 100)


#### pdu_M1 island ####
# Get tree path and abricate alignment data
path_to_tree <- "data/iqtree/core_gene_alignment.OGremoved.contree"
path_to_abricate <- "data/islandviewer/abricate_screening/ST372.pdu_M1.txt"
plasrefname <- "pdu_M1"

# Define reference length so we can calculate %ID later
pdu_M1_ref_len <- read_delim(path_to_abricate, n_max = 1) %>% 
  mutate(COVERAGE = as.integer(gsub("^.*\\/", "", COVERAGE))) %>% pull(COVERAGE)

# Run plasmidmapR to generate binned hit table for each ST372 sequence against pUTI89 sequence
plasmid_mapR(path_to_abricate = path_to_abricate,
             plasmid_reference_name = plasrefname,
             output_directory = NULL,
             min_hit_id = 90,
             min_hit_length = 0.5,
             writecsv = FALSE,
             path_to_tree = path_to_tree)

# Generate dataframe with a summary statistic for how many bins were covered at >90%ID
pdu_M1_ID <- as.data.frame(rowSums(pdu_M1_binned_hits))

# Add names column
pdu_M1_ID$Name <- rownames(pdu_M1_binned_hits)

# Rename columns
colnames(pdu_M1_ID) <- c("pdu_M1_ID","Name")

# Convert value to a percentage
pdu_M1_ID$pdu_M1_ID <- round((pdu_M1_ID$pdu_M1_ID/pdu_M1_ref_len) * 100)

#### kps_M2 island ####
# Get tree path and abricate alignment data
path_to_tree <- "data/iqtree/core_gene_alignment.OGremoved.contree"
path_to_abricate <- "data/islandviewer/abricate_screening/ST372.kps_M2.txt"
plasrefname <- "kps_M2"

# Define reference length so we can calculate %ID later
kps_M2_ref_len <- read_delim(path_to_abricate, n_max = 1) %>% 
  mutate(COVERAGE = as.integer(gsub("^.*\\/", "", COVERAGE))) %>% pull(COVERAGE)

# Run plasmidmapR to generate binned hit table for each ST372 sequence against pUTI89 sequence
plasmid_mapR(path_to_abricate = path_to_abricate,
             plasmid_reference_name = plasrefname,
             output_directory = NULL,
             min_hit_id = 90,
             min_hit_length = 0.5,
             writecsv = FALSE,
             path_to_tree = path_to_tree)

# Generate dataframe with a summary statistic for how many bins were covered at >90%ID
kps_M2_ID <- as.data.frame(rowSums(kps_M2_binned_hits))

# Add names column
kps_M2_ID$Name <- rownames(kps_M2_binned_hits)

# Rename columns
colnames(kps_M2_ID) <- c("kps_M2_ID","Name")

# Convert value to a percentage
kps_M2_ID$kps_M2_ID <- round((kps_M2_ID$kps_M2_ID/kps_M2_ref_len) * 100)

#### kps_M3 island ####
# Get tree path and abricate alignment data
path_to_tree <- "data/iqtree/core_gene_alignment.OGremoved.contree"
path_to_abricate <- "data/islandviewer/abricate_screening/ST372.kps_M3.txt"
plasrefname <- "kps_M3"

# Define reference length so we can calculate %ID later
kps_M3_ref_len <- read_delim(path_to_abricate, n_max = 1) %>% 
  mutate(COVERAGE = as.integer(gsub("^.*\\/", "", COVERAGE))) %>% pull(COVERAGE)

# Run plasmidmapR to generate binned hit table for each ST372 sequence against pUTI89 sequence
plasmid_mapR(path_to_abricate = path_to_abricate,
             plasmid_reference_name = plasrefname,
             output_directory = NULL,
             min_hit_id = 90,
             min_hit_length = 0.5,
             writecsv = FALSE,
             path_to_tree = path_to_tree)

# Generate dataframe with a summary statistic for how many bins were covered at >90%ID
kps_M3_ID <- as.data.frame(rowSums(kps_M3_binned_hits))

# Add names column
kps_M3_ID$Name <- rownames(kps_M3_binned_hits)

# Rename columns
colnames(kps_M3_ID) <- c("kps_M3_ID","Name")

# Convert value to a percentage
kps_M3_ID$kps_M3_ID <- round((kps_M3_ID$kps_M3_ID/kps_M3_ref_len) * 100)

####-----GENOMIC ISLAND FIGURES-----####
# Create a number of extra meta columns so they show up in the heatmap
bulk_meta <- meta %>% select(Name, Cluster)
bulk_meta <- cbind(bulk_meta, rep(bulk_meta[2], 19))
bulk_meta <- cbind(bulk_meta, meta %>% select(Source))
bulk_meta <- cbind(bulk_meta, rep(bulk_meta[22], 19))

#### pdu G island ####
# Convert values to binary categories
pdu_G_binned_hits[pdu_G_binned_hits < 100] <- "No"
pdu_G_binned_hits[pdu_G_binned_hits == 100] <- "Yes"

# In order to cbind the meta and the binned_hits we need to reorder the binned_hits rows as per the meta because they are not sorted the same way
pdu_G_binned_hits <- pdu_G_binned_hits[order(match(rownames(pdu_G_binned_hits), bulk_meta$Name)),]

# Combine hacked metadata and binned hits dataframe
pdu_G_binned_hits <- cbind(bulk_meta[,2:ncol(bulk_meta)], pdu_G_binned_hits)

# Plot the tree
gheatmap(scoary_tree, 
         data = pdu_G_binned_hits,
         font.size = 2,
         hjust = 0,
         colnames =FALSE,
         width = 20,
         offset = ,
         color = NULL) + 
  scale_fill_manual(values = c(bap_clrs, source_clrs, "Yes"="#EB8E8B", "No"="white", na.value="white")) +
  theme(legend.position = "none")

# Scale factor for plotting
scalefactor <- 2

# Copy plot 
dev.copy(png, res = 300*scalefactor, height=2480, width=3508, units = , 'outputs/figures/FigureS1.pduG.map.png')

# Save plot and dev off
dev.off()

#### pdu_M1 island ####
# Convert values to binary categories
pdu_M1_binned_hits[pdu_M1_binned_hits < 100] <- "No"
pdu_M1_binned_hits[pdu_M1_binned_hits == 100] <- "Yes"

# In order to cbind the meta and the binned_hits we need to reorder the binned_hits rows as per the meta because they are not sorted the same way
pdu_M1_binned_hits <- pdu_M1_binned_hits[order(match(rownames(pdu_M1_binned_hits), bulk_meta$Name)),]

# Combine hacked metadata and binned hits dataframe
pdu_M1_binned_hits <- cbind(bulk_meta[,2:ncol(bulk_meta)], pdu_M1_binned_hits)

# Plot the tree
gheatmap(scoary_tree, 
         data = pdu_M1_binned_hits,
         font.size = 2,
         hjust = 0,
         colnames =FALSE,
         width = 20,
         offset = ,
         color = NULL) + 
  scale_fill_manual(values = c(bap_clrs, source_clrs, "Yes"="#8DD3C7", "No"="white", na.value="white")) +
  theme(legend.position = "none")

# Copy plot 
dev.copy(png, res = 300*scalefactor, height=2480, width=3508, units = , 'outputs/figures/FigureS2.pduM1.map.png')

# Save plot and dev off
dev.off()


#### kps_M2 island ####
# Convert values to binary categories
kps_M2_binned_hits[kps_M2_binned_hits < 100] <- "No"
kps_M2_binned_hits[kps_M2_binned_hits == 100] <- "Yes"

# This one is a special case because a lot of sequences have no hits to it i.e. binned_hits is only 115 rows
# It can just be joined to meta and doesn't need the rowname rearrangement
kps_M2_binned_hits$Name <- rownames(kps_M2_binned_hits)

# Join to metadata
kps_M2_binned_hits <- left_join(meta %>%select(Name), kps_M2_binned_hits) %>% select(-Name)

# Convert NAs to 'No'
kps_M2_binned_hits[is.na(kps_M2_binned_hits)] <- "No"

# Combine hacked metadata and binned hits dataframe
kps_M2_binned_hits <- cbind(bulk_meta[,2:ncol(bulk_meta)], kps_M2_binned_hits)

# Give it rownames from meta (this is the correct order because it was joined to meta before)
rownames(kps_M2_binned_hits) <- meta$Name

# Plot the tree
gheatmap(scoary_tree, 
         data = kps_M2_binned_hits,
         font.size = 2,
         hjust = 0,
         colnames =FALSE,
         width = 20,
         offset = ,
         color = NULL) + 
  scale_fill_manual(values = c(bap_clrs, source_clrs, "Yes"="#8DD3C7", "No"="white"), na.value = "#000000") +
  theme(legend.position = "none")

# Copy plot 
dev.copy(png, res = 300*scalefactor, height=2480, width=3508, units = , 'outputs/figures/FigureS3.kpsM2.map.png')

# Save plot and dev off
dev.off()

#### kps_M3 island ####
# Convert values to binary categories
kps_M3_binned_hits[kps_M3_binned_hits < 100] <- "No"
kps_M3_binned_hits[kps_M3_binned_hits == 100] <- "Yes"

# This one is a special case because a lot of sequences have no hits to it i.e. binned_hits is only 115 rows
# It can just be joined to meta and doesn't need the rowname rearrangement
kps_M3_binned_hits$Name <- rownames(kps_M3_binned_hits)

# Join to metadata
kps_M3_binned_hits <- left_join(meta %>%select(Name), kps_M3_binned_hits) %>% select(-Name)

# Convert NAs to 'No'
kps_M3_binned_hits[is.na(kps_M3_binned_hits)] <- "No"

# Combine hacked metadata and binned hits dataframe
kps_M3_binned_hits <- cbind(bulk_meta[,2:ncol(bulk_meta)], kps_M3_binned_hits)

# Give it rownames from meta (this is the correct order because it was joined to meta before)
rownames(kps_M3_binned_hits) <- meta$Name

# Plot the tree
gheatmap(scoary_tree, 
         data = kps_M3_binned_hits,
         font.size = 2,
         hjust = 0,
         colnames =FALSE,
         width = 20,
         offset = ,
         color = NULL) + 
  scale_fill_manual(values = c(bap_clrs, source_clrs, "Yes"="#8DD3C7", "No"="white", na.value = "white")) +
  theme(legend.position = "none")

# Copy plot 
dev.copy(png, res = 300*scalefactor, height=2480, width=3508, units = , 'outputs/figures/FigureS4.kpsM3.map.png')

# Save plot and dev off
dev.off()


####-----HEATMAP PROCESSING-----####
# Split ABRicate gene hits into their functional groups
# Get all the hits from CARD database and intI1 and intI2. Fix all the messy names. Filter out genes present in >90% of isolates. These are housekeeping
# genes that sometimes mutate to confer AMR phenotypes but we are only concerned wiht acquired resistance genes
args <- geno_meta %>%
  select(all_of(meta_cols), starts_with("card"), contains("intI")) %>%
  rename_with(~ gsub("card_", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("Escherichia_coli_", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("_beta-lactamase", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("(", "_", .x, fixed = TRUE)) %>%
  rename_with(~ gsub(")", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("PC1__", "PC1_", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("EC_custom_intI1.*", "intI1", .x)) %>%
  rename_with(~ gsub("EC_custom_intI2.*", "intI2", .x)) %>%
  rename_with(~ gsub("Shigella_flexneri_chloramphenicol_acetyltransferase", "catA1", .x, fixed = TRUE)) %>% 
  select(where(is.character), where( ~ is.integer(.x) && sum(.x) <= .9*nrow(geno_meta))) %>%
  select(sort(names(.)), -ugd) %>%
  relocate(all_of(meta_cols), contains("intI"))

# Calculate total carriage of each gene in the collection
arg_totals <- t(args %>% summarise(across(where(is.integer), sum)))

arg_totals <- as_tibble(arg_totals, rownames = "Gene", .name_repair = "minimal") %>% 
  rename(Total = 2) %>% 
  mutate(Percentage = round(Total/nrow(geno_meta)*100, 2))

# Get hits from VFDB
vags <- geno_meta %>% 
  select(all_of(meta_cols), starts_with("vfdb")) %>%
  rename_with(~ gsub("vfdb_", "", .x))

# Select additional virulence genes from our custom database
custom_vags <- geno_meta %>% 
  select(Name, contains("EC_custom")) %>%
  rename_with(~ gsub("EC_custom_", "", .x, fixed = TRUE)) %>%
  mutate(usp = as.integer(rowSums(select(.,starts_with("usp"))))) %>%
  select(-starts_with("usp_")) %>%
  mutate(eitA = as.integer(rowSums(select(.,starts_with("eitA"))))) %>%
  select(-starts_with("eitA_")) %>% 
  select(Name,
         starts_with(c("cba", "cbi", "cjr", 
                       "cva", "cvi","eit",
                       "fecA", "hek", "hyl",
                       "iha","iss", "merA",
                       "ompT", "silA",
                       "terA", "traT", "usp"))) %>%
  rename_with(~ gsub("_[A-Z]{1,2}.*", "", .x)) %>%
  rename_with(~ gsub("_pAPEC-O1-ColBM", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_pUTI89", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_type3", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_|VFG1539", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_chromosomal", "_1", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_episomal", "_2", .x, fixed =TRUE))

# Join VFDB and custom hits and filter genes in less than 5% of strains
vags <- left_join(vags, custom_vags) %>% select(sort(names(.))) %>%
  relocate(all_of(meta_cols)) %>% 
  select(where(is.character), where( ~ is.integer(.x) && sum(.x) >= .05*nrow(vags)))

# Filter out genes we don't want due to rarity
# cps, fae, nle, yag
# vags_preheat <- vags %>% select(-starts_with(c("cps", "fae", "nle", "yag")))

# Split out operons we want to filter to marker genes
# clbA, entB, escC, espL*, espX*, espY*, fepA, fimH, gspM, shuA, tssA, ybtA
# Select the markers
vag_operon_markers <- vags %>% 
  select(Name, clbA, entB, starts_with("espY"), fepA, fimH, gspM,ybtA)

# Remove the operons to be replaced with markers
vags_preheat <- vags %>% 
 select(-starts_with(c("clb", "ent", "esc", "esp", "fep", "fim", "gsp", "shu", "tss", "ybt")))

# Rejoin the markers and filter out genes present in less than 10% of genomes
vags_preheat <- left_join(vags_preheat, vag_operon_markers) %>% 
   select(where(is.character), where( ~ is.integer(.x) && sum(.x) >= .10*nrow(vags))) %>% 
   select(sort(names(.))) %>%
   relocate(all_of(meta_cols))

# Show differences in VAGs when filtering out columns by total presence
# vag_cols <- names(vags)
# vag_cols5 <- names(vags %>% select(where(is.character), where( ~ is.integer(.x) && sum(.x) >= .05*nrow(vags))))
# vag_cols10 <- names(vags %>% select(where(is.character), where( ~ is.integer(.x) && sum(.x) >= .1*nrow(vags))))
# 
# setdiff(vag_cols, vag_cols5)
# setdiff(vag_cols5, vag_cols10)

# Calculate total carriage of each gene in the collection
vag_totals <- t(vags %>% summarise(across(where(is.integer), sum)))
vag_totals <- as_tibble(vag_totals, rownames = "Gene", .name_repair = "minimal") %>% 
  rename(Total = 2) %>% 
  mutate(Percentage = round(Total/nrow(geno_meta)*100, 2))

# Extract mobile genetic elements from the ISFinder database. NB: this data is not covered in the paper
mges <- geno_meta %>% 
  select(all_of(meta_cols), starts_with("ISfinder_")) %>%
  rename_with(~ gsub("ISfinder_Feb_2020_", "", .x)) %>%
  rename_with(~ gsub(":.*", "", .x))

# Extract plasmid related genes from the plasmidfinder database and tidy up messy names
plas <- geno_meta %>% 
  select(all_of(meta_cols), starts_with("plasmidfinder")) %>%
  rename_with(~ gsub("plasmidfinder_", "", .x)) %>%
  mutate(IncBOKZ = as.integer(rowSums(select(.,starts_with("IncB"))))) %>%
  select(-starts_with("IncB/")) %>%
  mutate(IncX1 = as.integer(rowSums(select(.,starts_with("IncX1"))))) %>%
  select(-starts_with("IncX1_")) %>%
  mutate(IncFII = as.integer(rowSums(select(.,starts_with("IncFII"))))) %>%
  select(-starts_with("IncFII("), -starts_with("IncFII_")) %>%
  mutate(IncFIA = as.integer(rowSums(select(.,starts_with("IncFIA"))))) %>%
  select(-starts_with("IncFIA("), -starts_with("IncFIA_")) %>%
  mutate(IncFIB = as.integer(rowSums(select(.,starts_with("IncFIB"))))) %>%
  select(-starts_with("IncFIB(")) %>%
  rename_with(~ gsub("_[0-9].*$", "", .x)) %>%
  rename_with(~ gsub("(", "_", .x, fixed = TRUE)) %>%
  rename_with(~ gsub(")", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("/", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("_FII", "", .x, fixed = TRUE)) %>%
  select(sort(names(.))) %>%
  relocate(all_of(meta_cols))

# Calculate total carriage of each gene in the collection
plas_totals <- t(plas %>% summarise(across(where(is.integer), sum)))
plas_totals <- as_tibble(plas_totals, rownames = "Gene") %>% 
  rename(Total = 2) %>% 
  mutate(Percentage = round(Total/nrow(geno_meta)*100, 2))

# Summary of IncF RSTs in the collection
# F_RST_totals <- geno_meta %>% 
#   filter(`F Plasmid` != "No F Plasmid") %>% 
#   group_by(`F Plasmid`, ColV, pUTI89) %>%
#   summarise(Total = n()) %>% 
#   mutate(Percentage = round(Total/nrow(geno_meta)*100, 2))

# Put dataframes into a list so they can be converted to matrices for heatmap visualisation
gene_list <- list(args = args, vags = vags_preheat, mges = mges, plas = plas)

# Loop that converts the gene screening dataframes into binary heatmaps named `geneprefix_heat` e.g args heatmap is called args_heat
for (f in 1:length(gene_list)){
  heat <- gene_list[[f]]
  heat <- heat %>% select(where(is_integer)) %>% as.matrix()
  rownames(heat) <- gene_list[[f]]$Name
  heat[is.na(heat)] <- 0
  heat[heat >= 1] <- "Present"
  heat[heat == 0] <- "Absent"
  heat <- as.data.frame(cbind(heat, meta %>% select(Source, Cluster)))
  heat <- heat %>% select(Source, Cluster, everything())
  assign(paste(names(gene_list[f]), "heat", sep="_"), heat)
}

# Additional colours for visualisation
heat_clrs <- c("Present" = "#000000E6", "Absent" = "#ededed", "Yes" = "#df03fc","No" = "white")

####-----CLUSTERS VS ARG/VAGs-----####
# Total ARGs
args2 <- args %>% select(-intI1, -intI2) %>% 
  mutate(`Total ARGs` = rowSums(across(where(is.integer)))) %>% 
  select(Name, Source, Cluster, `Total ARGs`)

# Average and maximum ARGs by Cluster
args2 %>% group_by(Cluster) %>% summarise(`Average ARGs` = mean(`Total ARGs`)) %>% arrange(desc(`Average ARGs`))
args2 %>% group_by(Cluster) %>% summarise(`Max ARGs` = max(`Total ARGs`)) %>% arrange(desc(`Max ARGs`))

# Average ARGs by Source
args2 %>% group_by(Source) %>% summarise(`Average ARGs` = mean(`Total ARGs`)) %>% arrange(desc(`Average ARGs`))

# Total VAGs
vags2 <- vags %>% mutate(`Total VAGs` = rowSums(across(where(is.integer)))) %>% 
  select(Name, Source, Cluster, `Total VAGs`)

# Average VAGs by Cluster
vags2 %>% group_by(Cluster) %>% summarise(`Average VAGs` = mean(`Total VAGs`)) %>% arrange(desc(`Average VAGs`))

# Average VAGs by Source
vags2 %>% group_by(Source) %>% summarise(`Average VAGs` = mean(`Total VAGs`)) %>% arrange(desc(`Average VAGs`))


####-----FIGURE S1-3 GENE HEATMAPS-----####
# Alignment of gene presence/absence with tree and metadata
# ggtree object for heatmap visualisation
gene_heatmap_tree <- ggtree(tree, branch.length = "none") %<+% meta

## AMR HEATMAP ##
figS5 <- gheatmap(gene_heatmap_tree, 
                     args_heat, 
                     width = 30,
                     font.size = 2,
                     colnames_offset_x = 0.00001,
                     colnames_offset_y = 2.2,
                     colnames_position = "top",
                     colnames_angle = 45,
                     hjust = 0,
                     color = NA) +
  scale_fill_manual(name = "Data", values = c(heat_clrs, bap_clrs, source_clrs)) +
  ylim(c(0,420)) +
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

## VIRULENCE HEATMAP ##
figS6 <- gheatmap(gene_heatmap_tree, 
                     vags_heat, 
                     width = 30,
                     font.size = 2,
                     colnames_offset_x = 0.00001,
                     colnames_offset_y = 2.2,
                     colnames_position = "top",
                     colnames_angle = 45,
                     hjust = 0,
                     color = NA) +
  scale_fill_manual(name = "Data", values = c(heat_clrs, bap_clrs, source_clrs)) +
  ylim(c(0,420)) +
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

## PLASMID REPLICON HEATMAP ##
figS7 <- gheatmap(gene_heatmap_tree, 
                     plas_heat, 
                     width = 30,
                     font.size = 2,
                     colnames_offset_x = 0.00001,
                     colnames_offset_y = 2.2,
                     colnames_position = "top",
                     colnames_angle = 45,
                     hjust = 0,
                     color = NA) +
  scale_fill_manual(name = "Data", values = c(heat_clrs, bap_clrs, source_clrs)) +
  ylim(c(0,420)) +
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

# Save Figure S5
ggsave("FigureS5_ARG_heatmap.pdf",
       figS5, 
       path = "outputs/figures/", 
       device = "pdf", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 300)

# Save Figure S6
ggsave("FigureS6_VAG_heatmap.pdf",
       figS6, 
       path = "outputs/figures/", 
       device = "pdf", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 300)

# Save Figure S7
ggsave("FigureS7_plas_heatmap.pdf",
       figS7, 
       path = "outputs/figures/", 
       device = "pdf", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = 300)

####-----TABLES S1, S2, S3-----####
## Process gene screening dataframes to generate Table S1
S1_args <- args %>% select(Name, where(is.integer))
S1_vags <- vags %>% select(Name, where(is.integer))
S1_plas <- plas %>% select(Name, where(is.integer))
S1_mges <- mges %>% select(Name, where(is.integer))

# Join them together
tableS1 <- left_join(S1_args, S1_vags)
tableS1 <- left_join(tableS1, S1_plas)
tableS1 <- left_join(tableS1, S1_mges)
tableS1 <- left_join(meta %>% rename(fimH_allele =fimH), tableS1)

# Write to file
write_csv(tableS1, "outputs/data/TableS1_fulldata.csv")

## Table S2 genes over-represented in fastbaps clusters
tableS2 <- rbind(cluster_M_over, cluster_G_over)
tableS2 <- rbind(tableS2, cluster_L_over)
tableS2 <- rbind(tableS2, cluster_J_over)

# Write to file
write_csv(tableS2, "outputs/data/TableS2_scoary.csv")

## Table S3 - BLAST alignment of MVC121 contig0001 to pdu operon from Salmonella (GenBank accession: AF026270.2)
tableS3 <- read_csv("data/blast/MVC121pdu_vs_Senterica_pdu.csv")

write_csv(tableS3, "outputs/data/TableS3_pduBLAST.csv")
