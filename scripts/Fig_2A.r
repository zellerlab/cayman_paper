
library(here)
library(readxl)
source(here('scripts', 'utils.r'))

# Load tax tree and define colors etc
tree.filtered <- read.tree(here('data', 'tree.genus.ncbi.filtered.nwk'))
taxa <- tree.filtered$node.label
unique_phyla <- taxa[str_detect(taxa, "p__")]
colors <- colorRampPalette(brewer.pal(max(12, length(unique_phyla)+2), 
                                      "Set2"))(length(unique_phyla)+2)
set.seed(3)
colors <- sample(colors, size = length(unique_phyla))

color_data <- data.frame(phyla = unique_phyla,
                         colors = colors)
color_data$colors[color_data$phyla == "p__Actinobacteria"] <- "#8286ed"
color_data$colors[color_data$phyla == "p__Tenericutes"] <- "#a65057"

substrate_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(8)[c(1,2,3,4,6)]


# Load tax profiles
# Load genomic CAZyme information on Almeida resource
almeidaCAZy <- read_tsv(here('data', 'almeida_cazy_annotations.tsv'))

# Load cazy annotations
#cazyAnnots <- read_tsv(here('data', "20230609_glycan_annotations_cleaned_Nic.tsv"))
completed_substrate_annotations <- read_xlsx(here("data", "Glycan_Annotations", "20230607_glycan_annotations_cleaned.xlsx"))
glycan_annotations_final_cleaned <- completed_substrate_annotations %>% select(c(Family:Subfamily,ORIGIN:FUNCTION_AT_DESTINATION_3, Glycan_annotation))
cazyAnnots <- glycan_annotations_final_cleaned

# Aggregate info for tree
almeidaCAZyAggregated <- almeidaCAZy %>%
  mutate(familyGroup = str_replace_all(almeidaCAZy$cazy_family, "[0-9_]", "")) %>%
  group_by(genome, Kingdom, Phylum, Class, Order, Family, Genus, Species, mOTU_ID, familyGroup, Genome_type) %>%
  summarize(numberCAZyGenes = length(unique(sequenceID)),
            numberUniqueCAZyFamilies = length(unique(cazy_family))) %>%
            #CAZyModuleDensityOverGenome = sum(annotationLength)/Length[1]) %>% 
  #pivot_longer(-c(genome, kingdom, class, class)) %>%
  pivot_longer(c(numberCAZyGenes, 
                 numberUniqueCAZyFamilies)) %>%
  rename(metric = name) %>%
  rename(species = Species,
      genus = Genus,
      family = Family,
      order = Order,
      class = Class,
      phylum = Phylum,
      kingdom = Kingdom) %>%
    #mutate(genus = map_chr(genus, \(x) str_c(str_split(x, " ")[[1]][2:length(str_split(x, " ")[[1]])], sep = " ", collapse = " "))) %>%
    mutate(genus = str_replace(genus, "^[a-zA-Z0-9]+ ", "")) %>%
  select(
    mOTU_ID,
    genus,
    familyGroup,
    metric,
    value
  ) %>%
  group_by(mOTU_ID, familyGroup, genus, metric) %>%
  summarize(value = mean(value)) %>%
  pivot_wider(id_cols = c(mOTU_ID, familyGroup, genus), 
              names_from = metric, 
              values_from = value) %>% 
  group_by(genus, familyGroup) %>% 
  summarize(N = length(numberCAZyGenes), 
            numberCAZyGenesSD = sd(numberCAZyGenes),
            numberCAZyGenes = mean(numberCAZyGenes), 
            numberUniqueCAZyFamiliesSD = sd(numberUniqueCAZyFamilies),            
            numberUniqueCAZyFamilies = mean(numberUniqueCAZyFamilies)) %>% 
  relocate(genus, 
           familyGroup,
           numberCAZyGenes, 
           numberCAZyGenesSD, 
           numberUniqueCAZyFamilies, 
           numberUniqueCAZyFamiliesSD,
           N) %>%
  filter(!is.na(genus)) %>%
  replace(is.na(.), 0) %>%
  mutate(genus = str_replace(genus, "g__", "")) %>%
  filter(familyGroup != "AA")
almeidaCAZyAggregated$familyGroup <- factor(almeidaCAZyAggregated$familyGroup, levels = cazyFamilyOrder)


file.copy(from = here('data', 'itol_template_files/dataset_symbols_template.txt'),
  here('data', 'itol_files/phylum_symbols_variation.txt'), overwrite = T)
filePath <- here('data', 'itol_files/phylum_symbols_variation.txt')
tx  <- readLines(filePath)
tx2  <- gsub(pattern = "MAXIMUM_SIZE,50", replace = "MAXIMUM_SIZE,25", x = tx)
writeLines(tx, con=filePath)

for (i in 1:dim(color_data)[1]){
  phylum <- color_data$phyla[i]
  write(str_c(phylum, "2", 10, "#000000", "1", "1",  sep = ",", collapse = ","), filePath, append = T)
}


# 2 Multi-Bar-chart indicating the family composition
file.copy(from = here('data', 'itol_template_files/dataset_multibar_template.txt'),
          here('data', 'itol_files/numberCAZYmes_multibar_variation.txt'), overwrite = T)
filePath <- here('data', 'itol_files/numberCAZYmes_multibar_variation.txt')
tx  <- readLines(filePath)
#tx2  <- gsub(pattern = "abc", replace = "ccccccccccccccccccccc", x = tx)
tx <- str_replace(tx, "#DATASET_SCALE,2000,10000,20000", "DATASET_SCALE,0-0-#000000-2-1-3, DATASET_SCALE,130-130-#000000-2-1-3,260-260-#000000-2-1-3,390-390-#000000-2-1-3")
tx <- str_replace(tx, "DATASET_LABEL,example multi bar chart", "DATASET_LABEL, Number of CAZymes found")
#tx <- str_replace(tx, "COLOR,#ff0000", "COLOR,#000000")
tx <- str_replace(tx, "FIELD_COLORS,#ff0000,#00ff00,#0000ff", str_c("FIELD_COLORS", str_c(familyColors, sep = ",", collapse = ','), collapse = ",", sep = ','))
tx <- str_replace(tx, "FIELD_LABELS,f1,f2,f3", str_c("FIELD_LABELS", str_c(levels(almeidaCAZyAggregated$familyGroup), sep = ",", collapse = ','), sep = ",", collapse = ","))
#tx <- str_replace(tx, "#WIDTH,1000", "WIDTH,163.5") # with manually setting 'Left margin' to -150, it will be perfectly aligned with the barplots
#tx <- str_replace(tx, "#BAR_SHIFT,0", "BAR_SHIFT,-150") # ... or perhaps programatically with this?
#tx <- str_replace(tx, "COLOR,#ff0000", "COLOR,#")
writeLines(tx, con=filePath)

tmp <- almeidaCAZyAggregated %>%
  group_by(genus) %>%
  arrange(match(familyGroup, levels(.$familyGroup))) %>%
  pivot_wider(id_cols = genus,
              names_from = familyGroup,
              values_from = numberCAZyGenes, values_fill = 0)
tmp$stri <- apply(tmp, 1, function(x) str_c(x[2:length(x)], sep = ",", collapse = ","))
tmp$stri <- str_replace_all(tmp$stri, " ", "")
tmp  <- tmp %>%
  filter(!is.na(genus))

for (i in 1:dim(tmp)[1]){
  family <- tmp$genus[i]
  write(str_c(family, tmp$stri[i], sep = ",", collapse = ","), filePath, append = T)
}

# 3 Multi-Bar-chart indicating the family composition
file.copy(from = here('data', 'itol_template_files/dataset_multibar_template.txt'),
  here('data', 'itol_files/numberuniqueCAZYmes_multibar_variation.txt'), overwrite = T)
filePath <- here('data', 'itol_files/numberuniqueCAZYmes_multibar_variation.txt')
tx  <- readLines(filePath)
#tx2  <- gsub(pattern = "abc", replace = "ccccccccccccccccccccc", x = tx)
tx <- str_replace(tx, "#DATASET_SCALE,2000,10000,20000", "DATASET_SCALE,0-0-#000000-2-1-3, DATASET_SCALE,35-35-#000000-2-1-3,70-70-#000000-2-1-3,105-105-#000000-2-1-3")
tx <- str_replace(tx, "DATASET_LABEL,example multi bar chart", "DATASET_LABEL, Number of unique CAZyme families found")
#tx <- str_replace(tx, "COLOR,#ff0000", "COLOR,#000000")
tx <- str_replace(tx, "FIELD_COLORS,#ff0000,#00ff00,#0000ff", str_c("FIELD_COLORS", str_c(familyColors, sep = ",", collapse = ','), collapse = ",", sep = ','))
tx <- str_replace(tx, "FIELD_LABELS,f1,f2,f3", str_c("FIELD_LABELS", str_c(levels(almeidaCAZyAggregated$familyGroup), sep = ",", collapse = ','), sep = ",", collapse = ","))
#tx <- str_replace(tx, "#WIDTH,1000", "WIDTH,163.5") # with manually setting 'Left margin' to -150, it will be perfectly aligned with the barplots
#tx <- str_replace(tx, "#BAR_SHIFT,0", "BAR_SHIFT,-150") # ... or perhaps programatically with this?
#tx <- str_replace(tx, "COLOR,#ff0000", "COLOR,#")
writeLines(tx, con=filePath)

tmp <- almeidaCAZyAggregated %>%
  group_by(genus) %>%
  arrange(match(familyGroup, levels(.$familyGroup))) %>%
  pivot_wider(id_cols = genus,
              names_from = familyGroup,
              values_from = numberUniqueCAZyFamilies, values_fill = 0)
tmp$stri <- apply(tmp, 1, function(x) str_c(x[2:length(x)], sep = ",", collapse = ","))
tmp$stri <- str_replace_all(tmp$stri, " ", "")
tmp  <- tmp %>%
  filter(!is.na(genus))

for (i in 1:dim(tmp)[1]){
  family <- tmp$genus[i]
  write(str_c(family, tmp$stri[i], sep = ",", collapse = ","), filePath, append = T)
}


#4 Heatmap with substrate z-scores
substrateNsByGenus <- almeidaCAZy %>% left_join(cazyAnnots %>%
                                                        select(Subfamily,
                                                               FUNCTION_AT_DESTINATION_1) %>%
                                                        rename(cazy_family = Subfamily) %>%
                                                        separate_rows(FUNCTION_AT_DESTINATION_1, sep = ','), by = 'cazy_family') %>%
  group_by(genome, mOTU_ID, Genus, FUNCTION_AT_DESTINATION_1) %>%
  tally()  %>%
  group_by(Genus, FUNCTION_AT_DESTINATION_1) %>%
  summarize(nSD = sd(n),
            n = mean(n)) %>%
  # NEW: remove Other/NA
  #filter(FUNCTION_AT_DESTINATION_1 != " Other" & !is.na(FUNCTION_AT_DESTINATION_1)) %>%
  filter(!is.na(FUNCTION_AT_DESTINATION_1)) %>%
  filter(!FUNCTION_AT_DESTINATION_1 == "Other") %>%
  filter(!is.na(Genus)) %>%
  replace(is.na(.), 0) %>% 
  pivot_wider(id_cols = Genus, names_from = FUNCTION_AT_DESTINATION_1, values_from = n, values_fill = 0)

# Subset to only what is there in tree before z-scoring (otherwise scales are getting all messed up...)
substrateNsByGenus$Genus <- str_replace(substrateNsByGenus$Genus, ".* ", "")
substrateNsByGenus <- substrateNsByGenus %>%
  filter(Genus %in% tree.filtered$tip.label)

# z-scores
substrateNsByGenus[, !colnames(substrateNsByGenus) == "Genus"] <- apply(substrateNsByGenus[, !colnames(substrateNsByGenus) == "Genus"], 2, scale)
# Cap values above 3 at 3/-3 at -3
substrateNsByGenus[, !colnames(substrateNsByGenus) == "Genus"] <- apply(substrateNsByGenus[, !colnames(substrateNsByGenus) == "Genus"], 2, function(x) ifelse(x > 3, 3, x))
substrateNsByGenus[, !colnames(substrateNsByGenus) == "Genus"] <- apply(substrateNsByGenus[, !colnames(substrateNsByGenus) == "Genus"], 2, function(x) ifelse(x < -3, -3, x))

substrateNsByGenus$stri <- apply(substrateNsByGenus, 1, function(x) str_c(as.character(x[2:length(x)]), sep = " ", collapse = " "))
substrateNsByGenus$stri <- str_replace_all(substrateNsByGenus$stri, "^ ", "")
substrateNsByGenus$stri <- str_replace_all(substrateNsByGenus$stri, "  ", " ")

file.copy(from = here('data', 'itol_template_files/dataset_heatmap_template.txt'),
          here('data', 'itol_files/dataset_heatmap_template.txt'), overwrite = T)
filePath <- here('data', 'itol_files/dataset_heatmap_template.txt')
tx  <- readLines(filePath)
tx  <- gsub(pattern = "#COLOR_MIN #ff0000", replace = "COLOR_MIN #0000FF", x = tx)
tx  <- gsub(pattern = "#COLOR_MAX #0000ff", replace = "COLOR_MAX #FF0000", x = tx)
tx  <- gsub(pattern = "#USE_MID_COLOR 1", replace = "USE_MID_COLOR 1", x = tx)
tx  <- gsub(pattern = "#COLOR_MID #ffff00", replace = "COLOR_MID #FFFFFF", x = tx)
tx  <- gsub(pattern = "FIELD_LABELS f1 f2 f3 f4 f5 f6", replace = str_c("FIELD_LABELS", str_c(substrateNsByGenus %>%
                                                                                                ungroup() %>%
                                                                                                 select(-c(Genus, stri)) %>%
                                                                                                 colnames(), sep = ' ', collapse = " "), sep = " ", collapse = " "), x = tx)
tx  <- gsub(pattern = "#USER_MIN_VALUE 0", replace = str_c("#USER_MIN_VALUE ", -6, sep = " ", collapse = " "), x = tx)
tx  <- gsub(pattern = "#USER_MID_VALUE 500", replace = "USER_MID_VALUE 0", x = tx)
tx  <- gsub(pattern = "#USER_MAX_VALUE 1000", replace = str_c("#USER_MAX_VALUE ", 6, sep = " ", collapse = " "), x = tx)
writeLines(tx, con=filePath)


for (i in 1:dim(substrateNsByGenus)[1]){
  family <- substrateNsByGenus$Genus[i]
  write(str_c(family, substrateNsByGenus$stri[i], sep = " ", collapse = " "), filePath, append = T)
}

# # 5 Leaf nodes scaled to the genus mean relative abundance
# file.copy(from = here('data', 'itol_template_files/dataset_symbols_template.txt'),
#           here('data', 'itol_files/dataset_symbols_template.txt'), overwrite = T)
# filePath <- here('data', 'itol_files/dataset_symbols_template.txt')
# tx  <- readLines(filePath)
# tx2  <- gsub(pattern = "MAXIMUM_SIZE,50", replace = "MAXIMUM_SIZE,25", x = tx)
# writeLines(tx, con=filePath)


# pPredGenus <- read_tsv(here('data', '2024_08_12_w_nw_profiles.tsv'))

# tmp <- pPredGenus %>% 
#   #select(-Unique_CAZyme_count) %>% 
#   # pivot_longer(-c(sampleID, 
#   #                 study_name, 
#   #                 non_westernized)) %>% 
#   #rename(genus = name) %>% 
#   ungroup() %>% 
#   group_by(genus) %>% 
#   #mutate(relAb = 10^value) %>% 
#   summarize(meanRelAb = mean(log10(relAb+pseudocount)))

# tmp  <- tmp %>%
#   filter(!is.na(genus))

# tmp <- tmp[tmp$genus %in% tree.filtered$tip.label, ]

# for (i in 1:dim(tmp)[1]){
#   family <- tmp$genus[i]
#   val <- tmp$meanRelAb[i[]]
#   write(str_c(family, "2", 1/-log10(val), "#000000", "1", "1",  sep = ",", collapse = ","), filePath, append = T)
# }

# Now upload tree.genus.ncbi.filtered.nwk + all generated itol_template_files to itol.embl.de to generate Fig 2A



