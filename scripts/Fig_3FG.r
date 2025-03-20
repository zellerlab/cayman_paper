library(here)
library(ComplexHeatmap)
library(patchwork)
library(ggplotify)
library(tidyverse)
library(readxl)
source(here('scripts', 'utils.r'))


phy <- c("Firmicutes",
"Bacteroidetes",
"Actinobacteria",
"Proteobacteria",
"Fusobacteria",
"Verrucomicrobia",
"Tenericutes",
"Spirochaetes",
"Euryarchaeota")

co <- c("#00BA38",
"#D39200",
"#F8766D",
"#00B9E3",
"#00C19F",
"#FF61C3",
"#DB72FB",
"#619CFF",
"#93AA00")

color_data_vec <- cbind(phy, co) %>%
  as.data.frame()
color_data_vec <- co
names(color_data_vec) <- phy


p <- read_tsv(here('data', '2024_08_12_w_nw_profiles.tsv'))

motus3.0_taxonomy <- read_tsv(here('data', "motus3.0_taxonomy.tsv"))

meta <- read_csv(here('data', '2024_08_12_Western_non_Western_metadata_meta_analysis.csv')) %>%
  rename(sampleID = sample_id)  

# Load cazy annotations
#cazyAnnots <- read_tsv(here('data', "20230609_glycan_annotations_cleaned_Nic.tsv"))
completed_substrate_annotations <- read_xlsx(here("data", "Glycan_Annotations", "20250219_Table_S1_incl_dbCAN3_annotations.xlsx"))
glycan_annotations_final_cleaned <- completed_substrate_annotations %>% select(c(Family:Subfamily,ORIGIN:FUNCTION_AT_DESTINATION_3, Glycan_annotation))
cazyAnnots <- glycan_annotations_final_cleaned

tree.filtered <- read.tree(here('data', "tree.genus.ncbi.filtered.nwk"))


#cazy_richness <- read_tsv(here('data', "20230929_CAZyme_count_df_Nic.tsv")) %>%
cazy_richness <- read_tsv(here('data', "Intermediate_Files", "20240809_CAZyme_count_df.tsv")) %>%  
  rename(sampleID = Sample) %>%
  select(all_of(c("sampleID", "Unique_CAZyme_count")))

almeidaCAZy <- read_tsv(here('data', 'almeida_cazy_annotations.tsv'))


#assocs_W <- read_csv(here('data', 'siamcat_NW_results_associations_for_Nic.csv')) %>%
assocs_W <- read_csv(here('data', "Intermediate_Files", 'siamcat_NW_results_associations.csv')) %>%  
  filter(p.adj < 0.05) %>%
  #rownames_to_column('cazyme_fam') %>%
  # select(famil,fc) %>%
  # mutate(topFamilyType = ifelse(fc < 0, "top_W_assoc", "top_NW_assoc")) %>%
  rename(family = cazyme_fam) %>%
  select(family, fc) %>%
  rename(fc_W = fc) %>%
  #filter(family %in% heatmapData$family) %>%
  arrange(desc(fc_W))

colAnnots <- cazyAnnots %>%
  select(Subfamily, FUNCTION_AT_DESTINATION_1) %>%
  distinct() %>%
  filter(!is.na(FUNCTION_AT_DESTINATION_1)) %>%
  separate_rows(FUNCTION_AT_DESTINATION_1, sep = ",") %>%
  mutate(m = 1) %>%
  group_by(Subfamily, FUNCTION_AT_DESTINATION_1) %>%
  distinct() %>%
  as.data.frame() %>%
  pivot_wider(id_cols = Subfamily, names_from = FUNCTION_AT_DESTINATION_1, values_from = m, values_fill = 0) %>%
  # Explicitly fill in NAs for families without any annotation
  full_join(data.frame(Subfamily = cazyAnnots$Subfamily) %>% distinct(), by = 'Subfamily') %>%
  #filter(Subfamily %in% colnames(basee)) %>%
  # left_join(cu, by = 'Subfamily') %>%
  # mutate(cazyFamilyCluster = as.factor(cazyFamilyCluster)) %>%
  # left_join(importantPairs, by = c("Subfamily" = "family")) %>%
  left_join(assocs_W, by = c("Subfamily" = "family")) %>%
  column_to_rownames("Subfamily") %>%
  rename(cazyme_enrichment = fc_W) %>%
  # mutate(cazyClass = map_chr(rownames(.), function(x) str_replace_all(x, '[0-9_]', ''))) %>%
  select(-c(Glycogen, Unknown, Other)) %>%
  filter(!is.na(cazyme_enrichment))

wEnrichedCAZymeFamilies <- colAnnots %>%
      filter(cazyme_enrichment < 0) %>%
      rownames(.)

fcs <- p %>%
    select(sampleID, taxon, study_name, non_westernized, relAb) %>%
    mutate(taxon = str_split_fixed(taxon, "[|]", n = 7)[, 6]) %>%
    mutate(taxon = str_replace(taxon, "g__", "")) %>%
    # Restrict analysis to genera shown in tree.
    inner_join(data.frame(taxon = tree.filtered$tip.label), by = 'taxon') %>%
    group_by(sampleID, taxon, study_name, non_westernized) %>%
    summarize(relAb = sum(relAb)) %>%
    pivot_wider(id_cols = c(sampleID, study_name, non_westernized), 
                names_from = taxon, 
                values_from = relAb, 
                values_fill = 0) %>%
    as.data.frame() %>%
    inner_join(cazy_richness, by = 'sampleID') %>%
    #select(-sampleID) %>%
    mutate(across(colnames(.)[!colnames(.) %in%  c('sampleID', 
                                                    'study_name', 
                                                    'Unique_CAZyme_count',
                                                    'non_westernized')], function(x) log10(x + 1E-5))) %>%
      select(-c(study_name, Unique_CAZyme_count)) %>%
          pivot_longer(-c(sampleID, non_westernized)) %>%
          as_tibble() %>%
          rename(genus = name, relAblog10 = value) %>%
          # mutate(relAb = 10^relAblog10 - 1E-5) %>%
          mutate(relAb = relAblog10) %>%
          group_by(genus) %>%
          nest() %>%
          mutate(gFC = map_dbl(data, CalcGFC)) %>%
          # right_join(data.frame(genus = hc1$labels)) %>%
          # mutate(genus = factor(genus, levels = rev(hc1$labels[hc1$order])))
          arrange(desc(gFC)) %>%
          # inner_join(data.frame(genus = unique(heatmapData$genus)))
          identity()
fcs$genus <- factor(fcs$genus, levels = fcs$genus)

uniqueFamiliesRaw <- almeidaCAZy %>%
  # We're only interested in Genera that we can name...
  filter(!str_detect(Genus, "gen[.]")) %>%
  filter(!str_detect(Genus, "incertae")) %>%
  group_by(Genus) %>%
  nest() %>%
  mutate(Genus = str_split_fixed(Genus, " ", n=2)[, 2]) %>% 
  mutate(presentFamilies = map(data, function(x) {
  numGenomes <- length(unique(x$genome))
  return(x %>%
    select(genome, cazy_family) %>%
    distinct() %>%
    group_by(cazy_family) %>%
    tally() %>%
    mutate(n_N = n / numGenomes) %>%
    # Call a family present if it's found in at least 20% of the genomes.
    filter(n_N > 0.20) %>%
    pull(cazy_family))
}))

data <- fcs %>%
  arrange(gFC) %>%
  select(-data) %>%
  mutate(gFC = -1 * gFC) %>%
  inner_join(
    uniqueFamiliesRaw %>%
      mutate(presentWenrichedFamilies = map(presentFamilies, \(x) x[x %in% wEnrichedCAZymeFamilies])) %>%
      mutate('W-enriched CAZymes genomically present' = map_dbl(presentWenrichedFamilies, length)) %>%
      rename(genus = Genus),
      by = 'genus'
  ) %>%
  left_join(motus3.0_taxonomy %>%
    select(Phylum, Genus) %>%
    mutate(phylum = str_split_fixed(Phylum, " ", n = 2)[, 2]) %>%
    mutate(genus = str_split_fixed(Genus, " ", n = 2)[, 2]) %>%
    select(genus, phylum) %>%
    distinct(), by = 'genus') %>%
    filter(gFC >= 0)
pLow <- ggplot() +
  geom_point(data = data, aes(x = gFC, y = `W-enriched CAZymes genomically present`, color = phylum)) +
  ggrepel::geom_text_repel(data = data %>% filter(`W-enriched CAZymes genomically present` > 75), aes(x = gFC, y = `W-enriched CAZymes genomically present`, label = genus, color = phylum)) +
  theme_classic() +
  #ggtitle('W-enriched CAZymes\ndata taken genomically') +
  scale_color_manual(values = color_data_vec) +
  ylab("western-enriched CAZymes") +
  xlab("enrichment in western [gFC]")


ggsave(plot = pLow, filename = here('figures', "Fig3_F.pdf"), width = 5, height = 2.85)  

###
###

prev_mOTUs <- p %>% 
  #group_by(sampleID) %>%
  #mutate(relAb = count / sum(count)) %>%
  ungroup() %>%
  group_by(taxon) %>%
  summarize(prevalence = mean(relAb > 0), maxAb = any(relAb > 0.01)) %>%
  filter(prevalence > 0.05) %>%
  filter(maxAb)  %>% 
  mutate(genus = str_split_fixed(taxon, "[|]", n= 7)[,  6])

pPred <- p %>% 
  select(sampleID, taxon, study_name, non_westernized, relAb) %>% 
  # for pPred, take only prevalent mOTUs
  inner_join(prev_mOTUs %>%
               select(taxon),
             by = 'taxon') %>%  
  pivot_wider(id_cols = c(sampleID, study_name, non_westernized), 
              names_from = taxon, 
              values_from = relAb, 
              values_fill = 0) %>%
  as.data.frame() %>%
    inner_join(cazy_richness, by = 'sampleID') %>%
    # select(-sampleID) %>%
    mutate(across(colnames(.)[!colnames(.) %in% c('sampleID',
      'study_name',
      'Unique_CAZyme_count',
      'non_westernized')], function(x) log10(x + 1E-5)))
colnames(pPred) <- map_chr(colnames(pPred), function(x) str_split(x, "[|]")[[1]][length(str_split(x, "[|]")[[1]])])  

# NOW WITH CORRECT ORDER
motus3.0_taxonomy_fixed <- motus3.0_taxonomy %>%
  mutate(mOTUs_ID = str_replace(mOTUs_ID, "_v3_", "_v31_")) %>%
  mutate(genus = str_split_fixed(Genus, " ", n = 2)[, 2])

pPredGenusPrevalenceCum <- pPred %>%
  select(-c(study_name, non_westernized, Unique_CAZyme_count)) %>%
  pivot_longer(-sampleID) %>%
  mutate(prevalent = value > -5) %>%
  rename(mOTUs_ID = name) %>%
  left_join(motus3.0_taxonomy %>% mutate(mOTUs_ID = str_replace(mOTUs_ID, "_v3_", "_v31_")) %>% select(mOTUs_ID, Genus), by = 'mOTUs_ID') %>%
  # losing 3
  filter(!is.na(Genus)) %>%
  group_by(sampleID, Genus) %>%
  summarize(cumulativePrevalence = sum(prevalent)) %>%
  mutate(genus = str_split_fixed(Genus, " ", n = 2)[, 2]) %>%
  select(-Genus) %>%
  left_join(pPred %>% select(sampleID, study_name, non_westernized, Unique_CAZyme_count), by = 'sampleID') %>%
  pivot_wider(id_cols = c(sampleID, study_name, non_westernized, Unique_CAZyme_count), names_from = genus, values_from = cumulativePrevalence) %>%
  ungroup()

pPredBase <- pPred
pPredGenusPrevalenceCum <- pPredGenusPrevalenceCum[match(pPredBase$sampleID, pPredGenusPrevalenceCum$sampleID), ]
stopifnot(all(pPredBase$sampleID == pPredGenusPrevalenceCum$sampleID))

currentGenera <- c()
allGenera <- c()
currentR2 <- c()
totalR2 <- c()
modelsTotal <- list()
modelsCurrent <- list()
numberSpecies <- list()
# generaLeft <- genera
#generaLeft <- colnames(pPredGenus)[! colnames(pPredGenus) %in% c("sampleID", "study_name", "non_westernized",'Unique_CAZyme_count')]
generaLeft <- tree.filtered$tip.label
generaLeft <- generaLeft[!generaLeft %in% c("Leuconostoc")]
i <- 1
while (TRUE){
  print("Selecting best genus group to add...")
  currentr2 <- 0
  for (ge in generaLeft){
    
    motus3.0_taxonomy_fixed_current <- motus3.0_taxonomy_fixed %>% filter(genus %in% ge) %>% pull(mOTUs_ID)
    pPredCurrent <- pPredBase %>% select(any_of(c('sampleID', "non_westernized", "Unique_CAZyme_count", motus3.0_taxonomy_fixed_current)))
    # Adding cum. mOTU prevalence as another feature makes practically no difference.
    pPredCurrent <- cbind(pPredCurrent, pPredGenusPrevalenceCum[ge])
    pPredCurrent <- pPredCurrent %>% select(-sampleID)
    currentModel <- lm(Unique_CAZyme_count ~ ., data = pPredCurrent %>% select(-non_westernized))
    if (summary(currentModel)$r.squared > currentr2) {
      currentr2 <- summary(currentModel)$r.squared
      currentGenus <- ge
      modelsCurrent[[i]] <- currentModel
      names(modelsCurrent)[i] <- ge
      currentR2[[i]] <- currentr2
      # Has to be changed due to the cum prevalence :)
      numberSpecies[[i]] <- dim(pPredCurrent)[2] - 3
    }
  }
  print(str_c("Adding ", currentGenus, " with individual R2 of ", currentr2))
  allGenera <- c(allGenera, currentGenus)
  blaaaTotal <- motus3.0_taxonomy_fixed %>% filter(genus %in% allGenera) %>% pull(mOTUs_ID)
  pPredTotal <- pPredBase %>% select(any_of(c("non_westernized", "Unique_CAZyme_count", blaaaTotal)))  
  totalModel <- lm(Unique_CAZyme_count ~ ., data = pPredTotal %>% select(-non_westernized))
  totalR2 <- c(totalR2, summary(totalModel)$r.squared)
  modelsTotal[[length(modelsTotal) + 1]] <- totalModel  

  i <- i + 1
  generaLeft <- generaLeft[!generaLeft == currentGenus]
  if (i == length(tree.filtered$tip.label)){
    break
  }  
}

data <- data.frame(genus = allGenera, R2cum = unlist(totalR2), R2indiv = unlist(currentR2), numSpecies = unlist(numberSpecies)) %>%
  as_tibble() %>%
  filter(R2indiv > 0.1) %>%
  left_join(motus3.0_taxonomy %>%
    select(Phylum, Genus) %>%
    mutate(phylum = str_split_fixed(Phylum, " ", n = 2)[, 2]) %>%
    mutate(genus = str_split_fixed(Genus, " ", n = 2)[, 2]) %>%
    select(-Phylum, -Genus) %>%
    distinct(), by = 'genus') %>%  
  mutate(genus = str_c(genus, " (", numSpecies, ")", sep = ""))

data$genus <- factor(data$genus, levels = data$genus)
p <- ggplot() +
  geom_line(data = data %>% pivot_longer(c(R2cum, R2indiv)) %>% filter(name == "R2cum"), aes(x = genus, y = value, group = 1)) +
  geom_bar(data = data %>% pivot_longer(c(R2cum, R2indiv)) %>% filter(name == "R2indiv"), aes(x = genus, y = value, fill = phylum), stat = 'identity') +
  #geom_text(data = data %>% pivot_longer(c(R2cum, R2indiv)) %>% filter(name == "R2indiv"), aes(x = genus, y = value - 0.02 * str_count(as.character(numSpecies)), label = numSpecies, angle = 90), stat = 'identity') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("R2") +
  scale_color_manual(color_data_vec) +
  NULL

#ggsave(plot = p, filename = "/g/scb2/zeller/karcher/CAZY_project_v2/plots/R2cumPlot_genus_wise_mOTU_abundances_added_KINDA_TRUE_FORWARD_SELECTION.pdf", width = 7, height = 3)
ggsave(plot = p, filename = here('figures', "Fig3_G.pdf"), width = 7, height = 3)