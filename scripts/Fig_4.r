library(here)
library(ComplexHeatmap)
library(patchwork)
library(ggplotify)
library(tidyverse)
library(readxl)
source(here('scripts', 'utils.r'))

# Running this whole script should take around half an hour to an hour.
# This script also produces extended data figure 7

cazyPc <- 0.01

p <- read_tsv(here('data', '2024_08_12_w_nw_profiles.tsv'))

# Load cazy annotations
#cazyAnnots <- read_tsv(here('data', "20230609_glycan_annotations_cleaned_Nic.tsv"))
completed_substrate_annotations <- read_xlsx(here("data", "Glycan_Annotations", "20230607_glycan_annotations_cleaned.xlsx"))
glycan_annotations_final_cleaned <- completed_substrate_annotations %>% select(c(Family:Subfamily,ORIGIN:FUNCTION_AT_DESTINATION_3, Glycan_annotation))
cazyAnnots <- glycan_annotations_final_cleaned

#cazy_richness <- read_tsv(here('data', "20230929_CAZyme_count_df_Nic.tsv")) %>%
cazy_richness <- read_tsv(here('data', 'Intermediate_Files', "20240809_CAZyme_count_df.tsv")) %>%  
  rename(sampleID = Sample) %>%
  select(all_of(c("sampleID", "Unique_CAZyme_count")))

motus3.0_taxonomy <- read_tsv(here('data', "motus3.0_taxonomy.tsv"))

meta <- read_csv(here('data', '2024_08_12_Western_non_Western_metadata_meta_analysis.csv')) %>%
  rename(sampleID = sample_id)  

almeidaCAZy <- read_tsv(here('data', 'almeida_cazy_annotations.tsv'))

tree.filtered <- read.tree(here('data', "tree.genus.ncbi.filtered.nwk"))

# Compute families present/genus
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

prev_mOTUs <- p %>% 
  #group_by(sampleID) %>%
  #mutate(relAb = count / sum(count)) %>%
  ungroup() %>%
  group_by(taxon) %>%
  summarize(prevalence = mean(relAb > 0), maxAb = any(relAb > 0.01)) %>%
  filter(prevalence > 0.05) %>%
  filter(maxAb)  %>% 
  mutate(genus = str_split_fixed(taxon, "[|]", n= 7)[,  6])

# Species-level profiles
pPredCountsWithoutPseudoCount <- p %>% 
  select(sampleID, taxon, study_name, non_westernized, count) %>% 
  # for pPred, take only prevalent mOTUs
  inner_join(prev_mOTUs %>%
               select(taxon),
             by = 'taxon') %>%  
  pivot_wider(id_cols = c(sampleID, study_name, non_westernized), 
              names_from = taxon, 
              #values_from = relAb, 
              values_from = count,
              values_fill = 0) %>%
  as.data.frame() %>%
  inner_join(cazy_richness, by = 'sampleID')

# Genus-level profiles
pPredGenusCountsWithoutPseudoCount <- p %>%
  select(sampleID, taxon, study_name, non_westernized, count) %>%
  mutate(taxon = str_split_fixed(taxon, "[|]", n = 7)[, 6]) %>%
  mutate(taxon = str_replace(taxon, "g__", "")) %>%
  # Restrict analysis to genera shown in tree.
  inner_join(data.frame(taxon = tree.filtered$tip.label), by = 'taxon') %>%
  group_by(sampleID, taxon, study_name, non_westernized) %>%
  summarize(count = sum(count)) %>%
  pivot_wider(id_cols = c(sampleID, study_name, non_westernized), 
              names_from = taxon, 
              values_from = count, 
              values_fill = 0) %>%
  as.data.frame() %>%
  inner_join(cazy_richness, by = 'sampleID')

colnames(pPredCountsWithoutPseudoCount) <- map_chr(colnames(pPredCountsWithoutPseudoCount), function(x) str_split(x, "[|]")[[1]][length(str_split(x, "[|]")[[1]])])  

pPredCountsWithoutPseudoCountCLR <- pPredCountsWithoutPseudoCount %>%
    relocate(sampleID, study_name, non_westernized, Unique_CAZyme_count)
pPredGenusCountsWithoutPseudoCountCLR <- pPredGenusCountsWithoutPseudoCount %>%
  relocate(sampleID, study_name, non_westernized, Unique_CAZyme_count)


pPredCountsWithoutPseudoCountCLR[, !colnames(pPredCountsWithoutPseudoCountCLR) %in% c("sampleID",
                                                                  "study_name",
                                                                  "non_westernized",
                                                                  "Unique_CAZyme_count")] <- 
  t(apply(pPredCountsWithoutPseudoCountCLR[, !colnames(pPredCountsWithoutPseudoCountCLR) %in% c("sampleID",
                                                                            "study_name",
                                                                            "non_westernized",
                                                                            "Unique_CAZyme_count")],
          # Do CLR (over samples!)
          # Also add literal pseudocount here
          1, function(x) compositions::clr(x + 1)@.Data))
pPredGenusCountsWithoutPseudoCountCLR[, !colnames(pPredGenusCountsWithoutPseudoCountCLR) %in% c("sampleID",
                                                                  "study_name",
                                                                  "non_westernized",
                                                                  "Unique_CAZyme_count")] <- 
  t(apply(pPredGenusCountsWithoutPseudoCountCLR[, !colnames(pPredGenusCountsWithoutPseudoCountCLR) %in% c("sampleID",
                                                                            "study_name",
                                                                            "non_westernized",
                                                                            "Unique_CAZyme_count")],
          # Do CLR (over samples!)
          # Also add literal pseudocount here
          1, function(x) compositions::clr(x + 1)@.Data))


#nwc <- readRDS(here('data', '20230929_complete_western_non_western_rel_abundances_wide_Nic.RDS'))
nwc <- readRDS(here('data', "Intermediate_Files", '20230929_complete_western_non_western_rel_abundances_wide.RDS'))

nwc <- nwc[, colnames(nwc) %in% pPredCountsWithoutPseudoCountCLR$sampleID]
pPredCountsWithoutPseudoCountCLR <- pPredCountsWithoutPseudoCountCLR[pPredCountsWithoutPseudoCountCLR$sampleID %in% colnames(nwc), ]
nwc <- nwc[, match(pPredCountsWithoutPseudoCountCLR$sampleID, colnames(nwc))]
#head(colnames(nwc) == pPredGenus$sampleID)
#stopifnot(all(colnames(nwc) == pPredGenus$sampleID))
pPredCountsWithoutPseudoCountCLR <- pPredCountsWithoutPseudoCountCLR[match(colnames(nwc), pPredCountsWithoutPseudoCountCLR$sampleID), ]
pPredGenusCountsWithoutPseudoCountCLR <- pPredGenusCountsWithoutPseudoCountCLR[match(colnames(nwc), pPredGenusCountsWithoutPseudoCountCLR$sampleID), ]
head(colnames(nwc) == pPredCountsWithoutPseudoCountCLR$sampleID)
stopifnot(all(colnames(nwc) == pPredCountsWithoutPseudoCountCLR$sampleID))
stopifnot(all(colnames(nwc) == pPredGenusCountsWithoutPseudoCountCLR$sampleID))

nwc <- nwc + cazyPc


# Remove families thar are not seen in more than 5% of samples
# Get presence/absence
nwcPrev <- (nwc > cazyPc) * 1
preve <- apply(nwcPrev, 1, mean)
preveGenera <- apply(nwcPrev, 2, mean)
#prevalentFamilies <- rownames(nwcPrev)[preve > 0.1]
nwc <- nwc[preve>0.20, ]

nwclog10 <- log10(nwc)

stopifnot(all(pPredCountsWithoutPseudoCountCLR$sampleID == colnames(nwclog10)))
stopifnot(all(pPredCountsWithoutPseudoCountCLR$sampleID == colnames(nwclog10)))
pPredCountsWithoutPseudoCountCLR_W <- pPredCountsWithoutPseudoCountCLR %>% filter(sampleID %in% 
(meta %>% filter(non_westernized == "Westernized") %>% pull(sampleID)))
pPredCountsWithoutPseudoCountCLR_NW <- pPredCountsWithoutPseudoCountCLR %>% filter(sampleID %in%
    (meta %>% filter(non_westernized == "Non_Westernized") %>% pull(sampleID)))

stopifnot(all(pPredGenusCountsWithoutPseudoCountCLR$sampleID == colnames(nwclog10)))
stopifnot(all(pPredGenusCountsWithoutPseudoCountCLR$sampleID == colnames(nwclog10)))
pPredGenusCountsWithoutPseudoCountCLR_W <- pPredGenusCountsWithoutPseudoCountCLR %>% filter(sampleID %in% 
(meta %>% filter(non_westernized == "Westernized") %>% pull(sampleID)))
pPredGenusCountsWithoutPseudoCountCLR_NW <- pPredGenusCountsWithoutPseudoCountCLR %>% filter(sampleID %in%
  (meta %>% filter(non_westernized == "Non_Westernized") %>% pull(sampleID)))

nwclog10_W <- nwclog10[, colnames(nwclog10) %in% (meta %>% filter(non_westernized == "Westernized") %>% pull(sampleID))]
nwclog10_NW <- nwclog10[, colnames(nwclog10) %in% (meta %>% filter(non_westernized == "Non_Westernized") %>% pull(sampleID))]

## 1
test2ClrOnCountsBySampleRPKMlog10_W <- get_pw_assocs(matrix1 =  pPredGenusCountsWithoutPseudoCountCLR_W, matrix2 = nwclog10_W, generaToTest =  tree.filtered$tip.label, cazyFamiliesToTest = rownames(nwclog10_NW))
test2ClrOnCountsBySampleRPKMlog10_W <- test2ClrOnCountsBySampleRPKMlog10_W %>%
  mutate(r_squared = map_dbl(ls, function(x) x$r.squared))
test2ClrOnCountsBySampleRPKMlog10_W$pval_taxon_adjusted <- p.adjust(test2ClrOnCountsBySampleRPKMlog10_W$pval, method = "BH")


## 2
test2ClrOnCountsBySampleRPKMlog10_NW <- get_pw_assocs(matrix1 =  pPredGenusCountsWithoutPseudoCountCLR_NW, matrix2 = nwclog10_NW, generaToTest =  tree.filtered$tip.label, cazyFamiliesToTest = rownames(nwclog10_NW))
test2ClrOnCountsBySampleRPKMlog10_NW <- test2ClrOnCountsBySampleRPKMlog10_NW %>%
  mutate(r_squared = map_dbl(ls, function(x) x$r.squared)) 
test2ClrOnCountsBySampleRPKMlog10_NW$pval_taxon_adjusted <- p.adjust(test2ClrOnCountsBySampleRPKMlog10_NW$pval, method = "BH")

## 3
liist <- list()
motus3.0_taxonomy$genus <- str_split_fixed(motus3.0_taxonomy$Genus, " ", n= 2)[,2]
for (g in tree.filtered$tip.label) {
   # get mOTU_IDs of genus
   m <- motus3.0_taxonomy %>% filter(genus == g) %>%
     mutate(mOTUs_ID = str_replace(mOTUs_ID, "_v3_", "_v31_")) %>%
     pull(mOTUs_ID)
  tmp <- pPredCountsWithoutPseudoCountCLR_W[colnames(pPredCountsWithoutPseudoCountCLR_W) %in% (m)]
  liist[[length(liist) + 1]] <- tmp
  names(liist)[length(liist)] <- g
}

testSpeciesFits_W <- get_pw_assocs_species_level_per_genus(liist, nwclog10_W, tree.filtered$tip.label, rownames(nwclog10_W))
testSpeciesFits_W <- testSpeciesFits_W %>%
  mutate(r_squared = map_dbl(ls, function(x) x$r.squared))
testSpeciesFits_W$pval_taxon_adjusted <- p.adjust(testSpeciesFits_W$pval, method = "BH")

# 4
liist <- list()
motus3.0_taxonomy$genus <- str_split_fixed(motus3.0_taxonomy$Genus, " ", n= 2)[,2]
for (g in tree.filtered$tip.label) {
   # get mOTU_IDs of genus
   m <- motus3.0_taxonomy %>% filter(genus == g) %>%
     mutate(mOTUs_ID = str_replace(mOTUs_ID, "_v3_", "_v31_")) %>%
     pull(mOTUs_ID)
  tmp <- pPredCountsWithoutPseudoCountCLR_NW[colnames(pPredCountsWithoutPseudoCountCLR_NW) %in% (m)]
  liist[[length(liist) + 1]] <- tmp
  names(liist)[length(liist)] <- g
}
testSpeciesFits_NW <- get_pw_assocs_species_level_per_genus(liist, nwclog10_NW, tree.filtered$tip.label, rownames(nwclog10_NW))
testSpeciesFits_NW <- testSpeciesFits_NW %>%
  mutate(r_squared = map_dbl(ls, function(x) x$r.squared))
testSpeciesFits_NW$pval_taxon_adjusted <- p.adjust(testSpeciesFits_NW$pval, method = "BH")

heatmapData_W_genus <- test2ClrOnCountsBySampleRPKMlog10_W %>%
  rename(beta_taxon = beta) %>%
  #filter(r_squared > 0.01) %>%
  #filter(beta_taxon >= 0) %>%
  #filter(pval_taxon_adjusted < 1E-5) %>%
  mutate(beta_taxon = ifelse(pval_taxon_adjusted < 0.01, beta_taxon, 0)) %>%
  inner_join(uniqueFamiliesRaw %>%
    select(Genus, presentFamilies) %>%
    unnest() %>%
    rename(genus = Genus, 
    family = presentFamilies), by = c('genus', 'family')) %>%
    left_join(motus3.0_taxonomy %>% 
    select(Phylum, Class, Order, Family, Genus) %>%
    filter(!apply(., 1, function(x) any(str_detect(x, "NA")))) %>%
    distinct() %>% 
    #mutate(across(all_of(colnames(.)), function(x) str_split_fixed(x, " ", n = 3)[, 2])) %>%
    mutate(Genus = str_split_fixed(Genus, " ", n = 3)[, 2]) %>%
    rename(genus = Genus), by = 'genus') %>%
    # select(genus, family, r_squared, beta_taxon, beta_W_offset, beta_W_interaction, pval_taxon_adjusted, pval_W_offset_adjusted, pval_W_interaction_adjusted, Phylum) %>%
    select(genus, family, r_squared, beta_taxon, pval_taxon_adjusted, Phylum)

heatmapData_NW_genus <- test2ClrOnCountsBySampleRPKMlog10_NW %>%
  rename(beta_taxon = beta) %>%
  #filter(r_squared > 0.01) %>%
  #filter(beta_taxon >= 0) %>%
  #filter(pval_taxon_adjusted < 1E-5) %>%
  mutate(beta_taxon = ifelse(pval_taxon_adjusted < 0.01, beta_taxon, 0)) %>%
  inner_join(uniqueFamiliesRaw %>%
    select(Genus, presentFamilies) %>%
    unnest() %>%
    rename(genus = Genus, 
    family = presentFamilies), by = c('genus', 'family')) %>%
    left_join(motus3.0_taxonomy %>% 
    select(Phylum, Class, Order, Family, Genus) %>%
    filter(!apply(., 1, function(x) any(str_detect(x, "NA")))) %>%
    distinct() %>% 
    #mutate(across(all_of(colnames(.)), function(x) str_split_fixed(x, " ", n = 3)[, 2])) %>%
    mutate(Genus = str_split_fixed(Genus, " ", n = 3)[, 2]) %>%
    rename(genus = Genus), by = 'genus') %>%
    # select(genus, family, r_squared, beta_taxon, beta_W_offset, beta_W_interaction, pval_taxon_adjusted, pval_W_offset_adjusted, pval_W_interaction_adjusted, Phylum) %>%
    select(genus, family, r_squared, beta_taxon, pval_taxon_adjusted, Phylum)

heatmapData_W <- testSpeciesFits_W %>%
#heatmapData_W <- test2ClrOnCountsBySampleRPKMlog10_W %>%
  rename(beta_taxon = beta) %>%
  #filter(r_squared > 0.01) %>%
  #filter(beta_taxon >= 0) %>%
  #filter(pval_taxon_adjusted < 1E-5) %>%
  mutate(beta_taxon = ifelse(pval_taxon_adjusted < 0.01, beta_taxon, 0)) %>%
  inner_join(uniqueFamiliesRaw %>%
    select(Genus, presentFamilies) %>%
    unnest() %>%
    rename(genus = Genus, 
    family = presentFamilies), by = c('genus', 'family')) %>%
    left_join(motus3.0_taxonomy %>% 
    select(Phylum, Class, Order, Family, Genus) %>%
    filter(!apply(., 1, function(x) any(str_detect(x, "NA")))) %>%
    distinct() %>% 
    #mutate(across(all_of(colnames(.)), function(x) str_split_fixed(x, " ", n = 3)[, 2])) %>%
    mutate(Genus = str_split_fixed(Genus, " ", n = 3)[, 2]) %>%
    rename(genus = Genus), by = 'genus') %>%
    # select(genus, family, r_squared, beta_taxon, beta_W_offset, beta_W_interaction, pval_taxon_adjusted, pval_W_offset_adjusted, pval_W_interaction_adjusted, Phylum) %>%
    select(genus, family, r_squared, beta_taxon, pval_taxon_adjusted, Phylum)
heatmapData_NW <- testSpeciesFits_NW %>%
#heatmapData_NW <- test2ClrOnCountsBySampleRPKMlog10_NW %>%
  rename(beta_taxon = beta) %>%
  #filter(r_squared > 0.01) %>%
  #filter(beta_taxon >= 0) %>%
  #filter(pval_taxon_adjusted < 1E-5) %>%
  mutate(beta_taxon = ifelse(pval_taxon_adjusted < 0.01, beta_taxon, 0)) %>%
  inner_join(uniqueFamiliesRaw %>%
    select(Genus, presentFamilies) %>%
    unnest() %>%
    rename(genus = Genus, 
    family = presentFamilies), by = c('genus', 'family')) %>%
    left_join(motus3.0_taxonomy %>% 
    select(Phylum, Class, Order, Family, Genus) %>%
    filter(!apply(., 1, function(x) any(str_detect(x, "NA")))) %>%
    distinct() %>% 
    #mutate(across(all_of(colnames(.)), function(x) str_split_fixed(x, " ", n = 3)[, 2])) %>%
    mutate(Genus = str_split_fixed(Genus, " ", n = 3)[, 2]) %>%
    rename(genus = Genus), by = 'genus') %>%
    # select(genus, family, r_squared, beta_taxon, beta_W_offset, beta_W_interaction, pval_taxon_adjusted, pval_W_offset_adjusted, pval_W_interaction_adjusted, Phylum) %>%
    select(genus, family, r_squared, beta_taxon, pval_taxon_adjusted, Phylum)


test <- inner_join(
    full_join(heatmapData_W_genus, heatmapData_NW_genus, by = c("genus", "family", "Phylum"), suffix = c(".westernized", ".non_westernized")),
    full_join(heatmapData_W, heatmapData_NW, by = c("genus", "family", "Phylum"), suffix = c(".westernized", ".non_westernized")), 
    by = c('genus', 'family'),
    suffix = c(".genus_level", ".species_level")
)

ij <- data.frame(genus = c(
    "Bifidobacterium",
    "Collinsella",
    "Bacteroides",
    "Prevotella",
    "Akkermansia"), family = c(
    "GH13_3",
    "GH13_30",
    "GH95",
    "GH95",
    "GT31"))

plo <- ggplot() +
    geom_point(data = test, aes(x = r_squared.westernized.genus_level, y = r_squared.westernized.species_level), color = 'red', alpha = 0.5) +
    geom_point(data = test %>% inner_join(ij), aes(x = r_squared.westernized.genus_level, y = r_squared.westernized.species_level), color = 'black', alpha = 1) +
    geom_point(data = test, aes(x = r_squared.non_westernized.genus_level, y = r_squared.non_westernized.species_level), color = 'blue', alpha = 0.5) +
    geom_point(data = test %>% inner_join(ij), aes(x = r_squared.non_westernized.genus_level, y = r_squared.non_westernized.species_level), color = 'black', alpha = 1) +
    geom_text_repel(data = test %>% inner_join(ij), aes(x = r_squared.westernized.genus_level, y = r_squared.westernized.species_level, label = str_c(genus, "__", family)), color = 'red') +
    geom_text_repel(data = test %>% inner_join(ij), aes(x = r_squared.non_westernized.genus_level, y = r_squared.non_westernized.species_level, label = str_c(genus, "__", family)), color = 'blue') +
    scale_size_continuous(limits = c(0, 1)) +
    theme_classic() +
    xlab("R2 (genus-level)") +
    ylab("R2 (species-level)")
ggsave(plot = plo, filename = here("figures", "Extended_Data_Fig7_r_squared_comparisons.pdf"), width = 6, height = 6)


heatmapData <- full_join(heatmapData_W, heatmapData_NW, by = c("genus", "family", "Phylum"), suffix = c(".westernized", ".non_westernized"))
heatmapData <- heatmapData %>%
  # Remove cells where both models are non-significant
  filter(pval_taxon_adjusted.westernized < 0.01 & pval_taxon_adjusted.non_westernized < 0.01) %>%
  ### Keep cells where both models return positive interaction
  #filter(beta_taxon.westernized >= 0 & beta_taxon.non_westernized >= 0) %>%
  # beta_taxon now corresponds to the mean fraction of coefficients > 0
  #filter(beta_taxon.westernized >= 0.25 & beta_taxon.non_westernized >= 0.25) %>%
  mutate(rsq_max = ifelse(r_squared.westernized > r_squared.non_westernized, "westernHigher", "nonWesternHigher")) %>%
  mutate(rsq_ratio = r_squared.westernized/r_squared.non_westernized) %>%
  mutate(r_squared = pmap_dbl(list(rsq_ratio, r_squared.westernized, r_squared.non_westernized), function(r, w, nw) {
    #if (r == "westernHigher") {
    if (r > 2) {      
      return(w)
    #} else if (r == "nonWesternHigher") {
    } else if (r < (1/2)) {      
      return(nw)
    } else {
      return(mean(c(w,nw)))
    }
  })) %>%
  mutate(cellGroup = case_when(
    rsq_ratio > 2 ~ "W_Higher",
    rsq_ratio < (1/2) ~ "NW_Higher",
    .default = "similar"))


a <- heatmapData %>%
    group_by(genus) %>%
    nest() %>%
    filter(map_lgl(data, function(x) {
      #return(any(x$r_squared > 0.25))
      return(any(x$r_squared.westernized > 0.4) | any(x$r_squared.non_westernized > 0.4))
    })) %>%
    select(genus)
b <- heatmapData %>%
    group_by(family) %>%
    nest() %>%
    filter(map_lgl(data, function(x) {
      #return(any(x$r_squared > 0.4))
      return(any(x$r_squared.westernized > 0.4) | any(x$r_squared.non_westernized > 0.4))
    })) %>%
    select(family)  

heatmapData <- heatmapData %>%
  inner_join(a) %>%
  inner_join(b)


#assocs_W <- read_csv(here('data', 'siamcat_NW_results_associations_for_Nic.csv')) %>%
assocs_W <- read_csv(here('data', "Intermediate_Files", 'siamcat_NW_results_associations.csv')) %>%
  filter(p.adj < 0.05) %>%
  #rownames_to_column('cazyme_fam') %>%
  # select(famil,fc) %>%
  # mutate(topFamilyType = ifelse(fc < 0, "top_W_assoc", "top_NW_assoc")) %>%
  rename(family = cazyme_fam) %>%
  select(family, fc) %>%
  rename(fc_W = fc) %>%
  filter(family %in% heatmapData$family) %>%
  arrange(desc(fc_W))
  # inner_join(heatmapData %>% select(family), by = 'family')

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
    inner_join(data.frame(genus = unique(heatmapData$genus)))
  fcs$genus <- factor(fcs$genus, levels = fcs$genus)

taxonFCBar <- fcs %>%
  ggplot() +
  # geom_pointrange(aes(
  #   x = genus,
  #   y = relAbMean,
  #   ymin = relAbMean - relAbSd,
  #   ymax = relAbMean + relAbSd,
  #   color = non_westernized),
  # position = position_dodge(width = 0.2)) +
  geom_bar(aes(
    x = genus,
    y = gFC
  ), stat = 'identity') +
  ylab("gFC (W/NW)") +
  coord_flip() +
    theme_classic()    

rsq_taxon <- heatmapData %>%
  pivot_wider(id_cols = c(genus), names_from = family, values_from = r_squared, values_fill = NA)  
cell_group <- heatmapData %>%
  pivot_wider(id_cols = c(genus), names_from = family, values_from = cellGroup, values_fill = NA)    

te <- map(list(rsq_taxon, cell_group), function(x){
  gf <- x$genus
  x <- x %>%
  as.data.frame() %>%
  select(-genus) %>%
  as.matrix()
  rownames(x) <- gf
  print(dim(x))
  x <- x[, colnames(x) %in% assocs_W$family]
  x <- x[rownames(x) %in% fcs$genus, ]
  print(dim(x))
  x <- x[, match(assocs_W$family, colnames(x))]
  x <- x[rev(match(fcs$genus, rownames(x))), ]
  #x <- x[,match(colnames(x), assocs_W$family)]
  return(x)
})
rsq_taxon <- te[[1]]
cell_group <- te[[2]]

basee <- rsq_taxon
tmp <- basee
tmp[is.na(tmp)] <- 0
hc1 <- hclust(dist(tmp, 'manhattan'))
hc2 <- hclust(dist(t(tmp), "manhattan"))


cu <- cutree(tree = hc2, k = 5) %>%
as.data.frame() %>%
rownames_to_column('Subfamily') %>%
rename(cazyFamilyCluster = ".")

myBreaks <- c(seq(min(basee, na.rm = T), 0, length.out=ceiling(100/2) + 1), 
              seq(max(basee, na.rm = T)/100, max(basee, na.rm = T), length.out=floor(100/2)))

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
  filter(Subfamily %in% colnames(basee)) %>%
  # left_join(cu, by = 'Subfamily') %>%
  # mutate(cazyFamilyCluster = as.factor(cazyFamilyCluster)) %>%
  # left_join(importantPairs, by = c("Subfamily" = "family")) %>%
  left_join(assocs_W, by = c("Subfamily" = "family")) %>%
  column_to_rownames("Subfamily") %>%
  rename(cazyme_enrichment = fc_W) %>%
  # mutate(cazyClass = map_chr(rownames(.), function(x) str_replace_all(x, '[0-9_]', ''))) %>%
  select(-c(Glycogen, Unknown, Other))
colAnnots <- colAnnots[match(colnames(basee), rownames(colAnnots)), ]

chosenOnes <- list(
    c("Bifidobacterium", "GH13_3"),
    c("Bifidobacterium", "GH13_30"),
    c("Collinsella", "GH13_30"),    
    c("Bacteroides", "GH95"),
    c("Prevotella", "GH95"),
    c("Akkermansia", "GT31")
)

chosenOnes <- map2(chosenOnes, str_to_upper(letters)[2:(length(chosenOnes) + 1)], function(x, num) {
  g <- x[1]
  f <- x[2]
  rowIndex <- which(rownames(basee) == g)
  colIndex <- which(colnames(basee) == f)
  if (any(is.null(list(g, f, rowIndex, colIndex)))) {
    print(str_c("Something wrong with ", g ,f, rowIndex, colIndex, sep = ", "))
  }
  return(c(x[1], x[2], rowIndex, colIndex, num))
})

taxaScatters <- map(chosenOnes, function(x) get_scatter(x[1], x[2], x[5]))
chosenOneIndices <- map(chosenOnes, function(x) c(x[3], x[4], x[5]))
chosenOneIndices <- data.frame(chosenOneIndices) %>% 
  t() %>% 
  as.data.frame() %>% 
  rename(x = V1, y = V2, ind = V3) %>% 
  mutate(across(all_of(c("x","y")), as.numeric))
l <- matrix(data = NA, nrow = dim(basee)[1], ncol = dim(basee)[2])
for (row in 1:dim(chosenOneIndices)[1]) {
  l[chosenOneIndices[row, 1], chosenOneIndices[row, 2]] <- chosenOneIndices[row, 3]
}
chosenOneIndices <- l

cell_fun_W = colorRamp2(c(0, max(basee[cell_group == "W_Higher" & !is.na(cell_group == "W_Higher")], na.rm = TRUE)), c("white", '#3B6FB6'))
cell_fun_NW = colorRamp2(c(0, max(basee[cell_group == "NW_Higher" & !is.na(cell_group == "NW_Higher")], na.rm = TRUE)), c("white", '#D41645'))
cell_fun_REST = colorRamp2(c(0, max(basee[cell_group == "similar" & !is.na(cell_group == "similar")], na.rm = TRUE)), c("white", 'purple'))

#pdf('/g/scb/zeller/karcher/tmp/test.pdf', width = 12, height = 10)
heatmapObjectLong <- Heatmap(
  as.matrix(t(basee)),
  # cluster_rows = hc1,
  cluster_rows = FALSE,
  # cluster_cols = hc2,
  cluster_columns = FALSE,
  #na_col = "lightgrey",
  na_col = "#CCCCCC",
  #na_col = "#E3E3E3", 
  column_names_gp = grid::gpar(fontsize = 12),
  row_names_gp = grid::gpar(fontsize = 6.5),
  column_split = factor(ifelse((fcs %>% .[dim(.)[1] : 1, ] %>% pull(gFC)) > 0, "NW-enriched\ngenera", "W-enriched\ngenera"), levels = c("W-enriched\ngenera", "NW-enriched\ngenera")),
  row_split = factor(ifelse(colAnnots$cazyme_enrichment > 0, "NW-enriched\ncazyme families", "W-enriched\ncazyme families"), levels = c("NW-enriched\ncazyme families", "W-enriched\ncazyme families")),
  #col = circlize::colorRamp2(c(min(basee, na.rm = TRUE), 0, max(basee, na.rm = TRUE)), c('#3B6FB6', "white", '#D41645')),
  col = circlize::colorRamp2(c(0, max(basee, na.rm = TRUE)), c("white", '#D41645')),  
  # left_annotation = rowAnnotation(`genus enrichment` = fcs[dim(fcs)[1] : 1, ]$gFC, col = list(`genus enrichment` = circlize::colorRamp2(c(min(fcs[dim(fcs)[1] : 1, ]$gFC), 0, max(fcs[dim(fcs)[1] : 1, ]$gFC)), c('#3B6FB6', "white", '#D41645')))),
  # top_annotation = rowAnnotation(`genus enrichment` = fcs[dim(fcs)[1] : 1, ]$gFC, col = list(`genus enrichment` = circlize::colorRamp2(c(min(colAnnots$cazyme_enrichment, na.rm = TRUE), 0, max(colAnnots$cazyme_enrichment, na.rm = TRUE)), c('#3B6FB6', "white", '#D41645')))),
  top_annotation = HeatmapAnnotation(`genus enrichment` = fcs[dim(fcs)[1]:1, ]$gFC,
    col = list(`genus enrichment` = circlize::colorRamp2(c(min(colAnnots$cazyme_enrichment, na.rm = TRUE), 0, max(colAnnots$cazyme_enrichment, na.rm = TRUE)),
      c('#3B6FB6', "white", '#D41645')))),
  left_annotation = rowAnnotation(
    `family enrichment` = colAnnots$cazyme_enrichment,
    DF = ifelse(is.na(colAnnots$DF), "NA", colAnnots$DF),
    PG = ifelse(is.na(colAnnots$PG), "NA", colAnnots$PG),
    GAG = ifelse(is.na(colAnnots$GAG), "NA", colAnnots$GAG),
    Mucin = ifelse(is.na(colAnnots$Mucin), "NA", colAnnots$Mucin),
    col = list(
      `family enrichment` = circlize::colorRamp2(c(min(colAnnots$cazyme_enrichment, na.rm = TRUE), 0, max(colAnnots$cazyme_enrichment, na.rm = TRUE)), c('#3B6FB6', "white", '#D41645')),
      `DF` = c("0" = "#CCCCCC", "1" = "black", "NA" = "white"),
      `PG` = c("0" = "#CCCCCC", "1" = "black", "NA" = "white"),
      `GAG` = c("0" = "#CCCCCC", "1" = "black", "NA" = "white"),
      `Mucin` = c("0" = "#CCCCCC", "1" = "black", "NA" = "white"))),
  # bottom_annotation = HeatmapAnnotation(
  #   W_full_model = W_full_model_fits$r_squared,
  #   NW_full_model = NW_full_model_fits$r_squared,
  #   r_square_diff = W_full_model_fits$r_squared - NW_full_model_fits$r_squared, 
  #   col = list(
  #     NW_full_model = circlize::colorRamp2(c(0, 1), c("white", '#D41645')),
  #     W_full_model = circlize::colorRamp2(c(0, 1), c("white", '#3B6FB6')),
  #     r_square_diff = circlize::colorRamp2(c(-0.5, 0, 0.5), c('#D41645', "white", '#3B6FB6'))
  #     )
  # ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (is.na(t(cell_group)[i, j])) {
        grid.rect(x, y, width, height, gp = gpar(fill = "#CCCCCC", col = "#CCCCCC"))
      } else if (t(cell_group)[i, j] == "W_Higher") {
        grid.rect(x, y, width, height, gp = gpar(fill = cell_fun_W(t(basee)[i, j]), col = cell_fun_W(t(basee)[i, j])))
      } else if (t(cell_group)[i, j] == "NW_Higher") {
        grid.rect(x, y, width, height, gp = gpar(fill = cell_fun_NW(t(basee)[i, j]), col = cell_fun_NW(t(basee)[i, j])))
      } else if (t(cell_group)[i, j] == "similar") {
        grid.rect(x, y, width, height, gp = gpar(fill = cell_fun_REST(t(basee)[i, j]), col = cell_fun_REST(t(basee)[i, j])))
        grid.rect(x, y, width, height, gp = gpar(fill = cell_fun_REST(t(basee)[i, j]), col = cell_fun_REST(t(basee)[i, j])))
      } else {
        asdasdasdasd
      }
      if (!is.na(t(chosenOneIndices)[i, j])){
        grid.text(t(chosenOneIndices)[i, j], x, y, gp = gpar(fontsize = 7, col = "white"))
      }
    }
  )
# heatmapObjectLong
# dev.off()

layout <- "
AAAAAAAAAAAAAAAHH
AAAAAAAAAAAAAAAHH
AAAAAAAAAAAAAAAHH
AAAAAAAAAAAAAAAHH
AAAAAAAAAAAAAAAHH
AAAAAAAAAAAAAAAHH
AAAAAAAAAAAAAAABB
AAAAAAAAAAAAAAABB
AAAAAAAAAAAAAAACC
AAAAAAAAAAAAAAACC
AAAAAAAAAAAAAAADD
AAAAAAAAAAAAAAADD
AAAAAAAAAAAAAAAEE
AAAAAAAAAAAAAAAEE
AAAAAAAAAAAAAAAFF
AAAAAAAAAAAAAAAFF
AAAAAAAAAAAAAAAGG
AAAAAAAAAAAAAAAGG
"


p <- wrap_plots(append(list(as.grob(heatmapObjectLong)), taxaScatters)) + patchwork::plot_layout(design = layout)

ggsave(plot = p, filename = here("figures", "Fig4.pdf"), width = 14, height = 15)
