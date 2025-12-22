library(tidyverse)
library(here)
library(ape)
library(ggembl)

source(here('scripts', 'utils.r'))
cazyPc <- 0.01

p <- read_tsv(here('data', '2024_08_12_w_nw_profiles.tsv'))

# total_bacterial_counts_per_sample <- p %>%
#   group_by(sampleID, study_name) %>%
#   summarize(total_bacterial_counts = sum(count))

motus3.0_taxonomy <- read_tsv(here('data', "motus3.0_taxonomy.tsv")) %>%
    rename(mOTU_ID = mOTUs_ID) %>%
    mutate(mOTU_ID = str_replace(mOTU_ID, "_v3_", "_v31_"))

meta <- read_csv(here('data', '2024_08_12_Western_non_Western_metadata_meta_analysis.csv')) %>%
  rename(sampleID = sample_id)  

tree.filtered <- read.tree(here('data', "tree.genus.ncbi.filtered.nwk"))

nwc <- readRDS(here('data', "Intermediate_Files", '20230929_complete_western_non_western_rel_abundances_wide.RDS'))
nwc <- nwc['GT31', ] %>% 
    as.data.frame() %>% 
    rename('GT31' = '.') %>%
    rownames_to_column('sampleID')

genus_level_profiles <- p %>%
    inner_join(motus3.0_taxonomy, by = 'mOTU_ID') %>%
  select(sampleID, taxon, study_name, non_westernized, count) %>%
  mutate(taxon = str_split_fixed(taxon, "[|]", n = 7)[, 6]) %>%
  mutate(taxon = str_replace(taxon, "g__", "")) %>%
  # Restrict analysis to genera shown in tree.
  inner_join(data.frame(taxon = tree.filtered$tip.label), by = 'taxon') %>%
  group_by(sampleID, taxon, study_name, non_westernized) %>%
  summarize(count = sum(count))

eukmot_files <- list.files('/g/scb/zeller/woerz/results/bio_studies/Cazy_Ducarmon/analysis/results_combined', full.names = T)
dataset_names <- str_remove(basename(eukmot_files), ".csv")
eukmot_profiles <- map2(eukmot_files, dataset_names, \(x, y) {
    read_csv(x) %>%
        mutate(dataset_name_patrick = y)
})

eukmot_profiles_long <- map(eukmot_profiles, \(x) {
    x %>%
        pivot_longer(-c(`species/samples`, dataset_name_patrick), names_to = "name", values_to = "value")
})

eukmot_profiles_long_df <- do.call(rbind, eukmot_profiles_long)
eukmot_profiles_long_df <- eukmot_profiles_long_df %>%
    filter(str_detect(`species/samples`, "Blasto")) %>%
    rename(taxon = `species/samples`, sampleID = name, count = value) %>%
    group_by(sampleID, dataset_name_patrick) %>%
    summarize(count = sum(count)) %>%
    rename(study_name = dataset_name_patrick) %>%
    ungroup()

#######################
#######################
#######################
#######################
#######################

library(tidyverse)
library(here)

sampleMetadata <- read_tsv('/g/scb/zeller/karcher/cayman_paper/data/cmd_meta.tsv')
cmd_meta <- sampleMetadata
cmd_meta_non_westernized <- cmd_meta %>% filter(body_site == "stool") %>% group_by(study_name) %>% filter(any(non_westernized == "yes"))
control_studies <- c("AsnicarF_2021", "HMP_2019_ibdmdb", "HMP_2019_t2d", "LeChatelierE_2013", "SchirmerM_2016", "ZeeviD_2015", "XieH_2016", "YachidaS_2019", "JieZ_2017", "QinJ_2012", "HMP_2012") 
profiled_controls_metadata <- cmd_meta %>% filter(body_site == "stool") %>% filter(study_name %in% control_studies) %>% filter(disease == "healthy")


########### 
# Liu FIX #
###########

Liu_2016_meta <- read_delim(here("data", "Metadata", "Liu_2016_sample_info"))
eukmot_profiles_long_df <- left_join(eukmot_profiles_long_df, Liu_2016_meta %>% select(run_accession, sample_alias), by = c("sampleID" = "sample_alias")) 
eukmot_profiles_long_df <- eukmot_profiles_long_df %>% mutate(run_accession = coalesce(run_accession, sampleID)) %>% select(-sampleID) %>% dplyr::rename("sampleID" = run_accession) ## dplyr::rename(run_accession = "sampleID") ## Maybe this should be dplyr::rename("sampleID" = run_accession)

################
# Schirmer FIX #
################
Schirmer_2016_metadata_SRA <- read_delim(here("data", "Metadata", "Schirmer_2016_SRA_Metadata.txt"))
## Remove all the single-end files from Schirmer, since these are all likely remnants of KNEADDATA pre-processing and are therefore supersmall files.
Schirmer_2016_metadata_SRA_single <- Schirmer_2016_metadata_SRA %>% filter(LibraryLayout == "SINGLE")
Schirmer_2016_metadata_SRA_single_vector <- Schirmer_2016_metadata_SRA_single$Run
## Filter out all Schirmer Single samples
#Western_df_list_rel_abun_long_merged <- Western_df_list_rel_abun_long_merged %>% filter(SampleID %nin% Schirmer_2016_metadata_SRA_single$Run)
#Schirmer_sample_names <- Western_df_list_rel_abun_long_merged %>% filter(Study_ID == "Schirmer_2016") %>% distinct(SampleID) %>% filter(SampleID != "unannotated")

profiled_controls_metadata_Schirmer <- profiled_controls_metadata %>% filter(study_name == "SchirmerM_2016")

## Filter out the double accession numbers (so throw out the singles) with the very nice function separate_rows from tidyverse
profiled_controls_metadata_Schirmer <- profiled_controls_metadata_Schirmer %>% print() %>% separate_rows(NCBI_accession, sep = ";") %>% ungroup() %>% select(
    NCBI_accession, sample_id
) %>% distinct()

eukmot_profiles_long_df <- left_join(eukmot_profiles_long_df, profiled_controls_metadata_Schirmer, by = c('sampleID' = "sample_id")) 
eukmot_profiles_long_df <- eukmot_profiles_long_df %>% mutate(run_accession = coalesce(NCBI_accession, sampleID)) %>% select(-sampleID) %>% dplyr::rename("sampleID" = run_accession) ## dplyr::rename(run_accession = "sampleID") ## Maybe this should be dplyr::rename("sampleID" = run_accession)

################
# Qin_2012 FIX #
################
# 1

eukmot_profiles_long_df <- eukmot_profiles_long_df %>%
    mutate(
        sampleID = ifelse(study_name == "qin_2012", str_replace(sampleID, "bgi-", ""), sampleID)
    )

genus_level_profiles %>% 
    group_by(
        sampleID, 
        study_name
        ) %>% 
        tally() %>% 
        left_join(
            eukmot_profiles_long_df %>% 
            select(-study_name) %>% 
            mutate(bla = 'bla'), 
        by = 'sampleID') %>%
        filter(is.na(bla)) %>% 
        group_by(study_name) %>% 
        tally()

sample_intersection <- intersect(
    genus_level_profiles$sampleID,
    eukmot_profiles_long_df$sampleID
); print(length(sample_intersection))

all_counts <- rbind(
    genus_level_profiles %>%
        ungroup() %>%
        inner_join(
            data.frame(sampleID = sample_intersection)
        ) %>%
        select(sampleID, study_name, taxon, count),
    eukmot_profiles_long_df %>%
        ungroup() %>%
        inner_join(
            data.frame(sampleID = sample_intersection)
        ) %>%    
        select(-NCBI_accession) %>%
        mutate(taxon = "Blastocystis") 
) %>%
    group_by(
        sampleID
    ) %>%
    mutate(relative_abundance = count / sum(count)) %>%
    mutate(total_counts_in_that_sample = sum(count)) %>%
    mutate(CLR_from_counts = compositions::clr(count + 1)@.Data) %>%
    mutate(CLR_from_rel_ab = compositions::clr(relative_abundance + (1/total_counts_in_that_sample))@.Data) %>%
    inner_join(
        data.frame(taxon = c("Blastocystis", "Akkermansia"))
    ) %>%
    select(
        sampleID,
        study_name,
        taxon,
        CLR_from_counts,
        CLR_from_rel_ab
    )    


all_data <- all_counts %>%
    left_join(
        meta %>% select(sampleID, non_westernized),  by = 'sampleID'
    ) %>%
    pivot_longer(
        -c(non_westernized, study_name, sampleID, taxon)
    ) %>%
    mutate(value  = ifelse(is.na(value), 0, value)) %>%
    #rename(taxon = name, taxon_abundance = value) %>%
    pivot_wider(
        id_cols = c(sampleID,taxon, study_name, non_westernized),
        names_from = name,
        values_from = c(value)
    ) %>%
    left_join(
        nwc %>% select(sampleID, GT31), by = 'sampleID'
    ) %>%
    mutate(GT31_log10 = log10(GT31 + cazyPc)) %>%
    select(-GT31)


pearson_cors <- all_data %>%
  mutate(`Income status` = ifelse(non_westernized == "Non_Westernized", "LMIS", "HIS")) %>%
  #group_by(taxon, non_westernized) %>%
    group_by(taxon, `Income status`) %>%
  summarize(
    n = sum(!is.na(CLR_from_rel_ab) & !is.na(GT31_log10)),
    cor = if (n >= 2) cor(CLR_from_rel_ab, GT31_log10, use = "complete.obs") else NA_real_,
    p_value = if (n >= 3) tryCatch(cor.test(CLR_from_rel_ab, GT31_log10)$p.value, error = function(e) NA_real_) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    r2 = cor^2,
    # formatted label with superscript 2 and two decimals (e.g. "R² = 0.45")
    label = ifelse(is.na(r2), "R² = NA", paste0("R\u00B2 = ", formatC(r2, format = "f", digits = 2)))
  )

p <- all_data %>%
  mutate(`Income status` = ifelse(non_westernized == "Non_Westernized", "LMIS", "HIS")) %>%
  ggplot() +
  geom_point(aes(x = CLR_from_rel_ab, y = GT31_log10, color = `Income status`), alpha = 0.2) +
  facet_grid(
    `Income status` ~ taxon,
    scales = 'free_x'
  ) +
  geom_text(
    data = pearson_cors,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1,
    vjust = 1.1,
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_embl() +
  scale_color_manual(
    name = "Income status",
    values = c("LMIS" = '#D41645', "HIS" = '#3B6FB6')
  ) +
  xlab("taxon abundance (CLR)") +
  ylab("GT31 (log10)") 

ggsave(
    plot = p,
    filename = here('figures', 'Akkermansia_Blastocystis_vs_GT31_abundance_CLR_log10.pdf'),
    width = 6,
    height = 4.5
)
