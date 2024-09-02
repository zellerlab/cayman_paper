library(tidyverse)
library(here)
library(readxl)
# library(pheatmap)
library(ComplexHeatmap)

pseudocount <- 1

completed_substrate_annotations <- read_xlsx(here("data", "Glycan_Annotations/", "20230607_glycan_annotations_cleaned.xlsx"))
glycan_annotations_final_cleaned <- completed_substrate_annotations %>% select(c(Family:Subfamily, ORIGIN:FUNCTION_AT_DESTINATION_3, Glycan_annotation))
cazyAnnots <- glycan_annotations_final_cleaned %>%
    separate_rows(FUNCTION_AT_DESTINATION_1, sep = ",") %>%
    select(Subfamily, FUNCTION_AT_DESTINATION_1) %>%
    mutate(present = ifelse(is.na(FUNCTION_AT_DESTINATION_1), 'NA', "Yes")) %>%
    pivot_wider(names_from = FUNCTION_AT_DESTINATION_1, values_from = present, values_fill = list(present = "No"))
cazyAnnots <- apply(cazyAnnots, 1, \(x) {
    if (x[2] != "No") {
        x[3:length(x)] <- "NA"
    }
    return(x)
})
cazyAnnots <- cazyAnnots %>% as.data.frame() %>% t() %>% as.data.frame() %>% select(-all_of(c("Other", "Unknown", "NA")))

source(here("scripts", "revisions_aditions", "explore_within_genus_cazy_diversity_utils.r"))

motus_level_agg <- read_tsv(here("data/motus_level_mean_cazy_abundances.tsv")) %>%
    # get mean numbers back
    mutate(across(all_of(colnames(.)[!get_base_col_name_boolean(colnames(.))]), \(x) str_split_fixed(x, " ", 2)[, 1])) %>%
    mutate(across(all_of(colnames(.)[!get_base_col_name_boolean(colnames(.))]), as.numeric)) %>%
    mutate(Species_simple = apply(str_split_fixed(Species, " ", n = 4)[, 2:3], 1, \(x) str_c(x, collapse = " "))) %>%
    mutate(mOTU_ID = str_split_fixed(mOTU_ID, "_", n = 5)[, 4]) %>%
    mutate(`log10_#_genomes` = log10(number_genomes)) %>%
    filter(number_genomes > 10)
genome_level <- read_tsv(here("data/genome_level_cazy_abundances.tsv")) %>%
    mutate(mOTU_ID = str_split_fixed(mOTU_ID, "_", n = 5)[, 4])


####################
# Compute genus-wise 'average variance per family' and look at most variable genera
# This is in direct response to the reviewer comment
# Therefore, showing which species contribute most CAZyme enriches would provide a more in-depth systematic overview that could be used as a resource by other researchers.
####################

# library(GGally)
# library(ggembl)
# dev_null <- motus_level_agg_metrics %>%
#     ungroup() %>%
#     select(
#         Genus,
#         # mean_variance_over_cazyme_z_scores,
#         mean_dist_to_other_motus,
#         variance_cazyme_z_scores,
#         sum_cazyme_copy_number,
#         sum_cazyme_families,
#         sum_unique_cazyme_families) %>%
#     unnest() %>%
#     group_by(Genus) %>%
#     nest() %>%
#     mutate(pairplots = map(data, \(x) {
#         x <- x %>% select(-all_of(c("Species", "mOTU_ID")))
#         return(ggpairs(x, upper = list(continuous = wrap(ggally_cor, method = "spearman"))))
#     })) %>%
#     mutate(genus_name_only = str_split_fixed(Genus, " ", 2)[, 2]) %>%
#     mutate(dev_null = pmap(list(genus_name_only, pairplots), \(genusname, p) {
#         genus_name_first_letter_uppercase <- str_to_title(genusname)
#         dir.create(here("figures", "revisions", str_c(genus_name_first_letter_uppercase, '_cazy_zoom_in')), showWarnings = FALSE)
#         ggsave(
#             filename = here("figures", "revisions", str_c(genus_name_first_letter_uppercase, '_cazy_zoom_in'), "metrics_pairplot.pdf"),
#             plot = p + theme_presentation(),
#             width = 10,
#             height = 10,
#             dpi = 300
#         )
#     }))


source(here("scripts", "revisions_aditions", "explore_within_genus_cazy_diversity_utils.r"))
get_whole_deal_motu_level(
    genus_number_and_name = "1263 Ruminococcus",
    taxa_to_add_explicitly = c("torques", "gnavus", "obeum"), data_transformation = "zscore"
)
get_whole_deal_motu_level(
    genus_number_and_name = "1263 Ruminococcus",
    taxa_to_add_explicitly = c("torques", "gnavus", "obeum"), data_transformation = "log10"
)
get_whole_deal_genome_level(
    taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "gnavus")) %>% pull(mOTU_ID),
    genus = "Ruminococcus",
    data_transformation = 'zscore',
    families_of_interest = NULL)
get_whole_deal_genome_level(
    taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "gnavus")) %>% pull(mOTU_ID),
    genus = "Ruminococcus",
    data_transformation = 'log10',
    families_of_interest = NULL)
get_whole_deal_genome_level(
    taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "torques")) %>% pull(mOTU_ID),
    genus = "Ruminococcus",
    data_transformation = 'zscore',
    families_of_interest = NULL)
get_whole_deal_genome_level(
    taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "torques")) %>% pull(mOTU_ID),
    genus = "Ruminococcus",
    data_transformation = 'log10',
    families_of_interest = NULL)
## We only have one species of Eisenbergiella, so we can't really do much here...
# get_whole_deal_motu_level(
#     genus_number_and_name = "1432051 Eisenbergiella",
#     taxa_to_add_explicitly = NULL
# )
source(here("scripts", "revisions_aditions", "explore_within_genus_cazy_diversity_utils.r"))
get_whole_deal_genome_level(
    taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "Eisenbergiella")) %>% pull(mOTU_ID),
    genus = "Eisenbergiella",
    data_transformation = 'zscore',
    families_of_interest = NULL)
get_whole_deal_genome_level(
    taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "Eisenbergiella")) %>% pull(mOTU_ID),
    genus = "Eisenbergiella",
    data_transformation = 'log10',
    families_of_interest = NULL)
## We only have one species of Hungatella with > 10 genomes, so we can't really do much here either...
# get_whole_deal_motu_level(
#     genus_number_and_name = "1432051 Eisenbergiella",
#     taxa_to_add_explicitly = NULL
# )
get_whole_deal_genome_level(
    taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "Hungatella")) %>% pull(mOTU_ID),
    genus = "Hungatella",
    data_transformation = 'zscore',
    families_of_interest = NULL)
get_whole_deal_genome_level(
    taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "Hungatella")) %>% pull(mOTU_ID),
    genus = "Hungatella",
    data_transformation = 'log10',
    families_of_interest = NULL)
get_whole_deal_motu_level(
    genus_number_and_name = "816 Bacteroides",
    data_transformation = 'zscore'
)
get_whole_deal_motu_level(
    genus_number_and_name = "816 Bacteroides",
    data_transformation = 'log10'
)
# This runs for too long - heatmap is very large
# get_whole_deal_genome_level(
#     taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "Bacteroides")) %>% filter(number_genomes > 500) %>% pull(mOTU_ID),
#     genus = "Bacteroides",
#     families_of_interest = NULL)
