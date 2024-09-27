##############
##############
### ATTENTION # Run from conda activate /g/scb/zeller/fspringe/Software/miniconda/envs/r_env_4.3.1
##############
##############

library(tidyverse)
library(here)
library(readxl)
# library(pheatmap)
library(ComplexHeatmap)
library(ggembl)
library(ggnewscale)
library(patchwork)
source(here('scripts', 'utils.r'))

pseudocount <- 1

# Mucin families shown in (old) Fig 2B, should now be Fig3A?
a <- c(
        "GH2",
        "GH92",
        "GH20",
        "GH31",
        "GH29",
        "GH95",
        "GH35",
        "GH33",
        "GH42",
        "GH130",
        "GH18",
        "GH109",
        "GH73",
        "GH84",
        "GH89",
        "GH123",
        "GH85",
        "GH112",
        "GH108"
        )

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

set.seed(1123)
family_variance_within_motus <- genome_level %>%
    group_by(mOTU_ID) %>%
    slice_sample(n = 50) %>%
    ungroup() %>%
    filter(!is.na(Genus)) %>%
    select(-c(
        Kingdom,
        Phylum,
        Class,
        Order
    )) %>%
    pivot_longer(
        -c(
            genome,
            Genome_type,
            Length,
            mOTU_ID,
            Family,
            Genus,
            Species
        ),
    ) %>%
    rename(
        cazy_family = name,
        copy_number = value
    ) %>%
    mutate(present = ifelse(copy_number > 0,  1, 0)) %>%
    group_by(mOTU_ID, Family, Genus, Species, cazy_family) %>%
    summarize(
        num_genomes = n(),
        prevalence = mean(present),
        mean_copy_num = mean(copy_number),
        # var = var(copy_number),
        # cov = sd(copy_number) / mean(copy_number)) %>%
        var = var(copy_number),
        cov = sd(present) / mean(present)) %>%        
    filter(num_genomes > 10)

only_mucin <- TRUE
all_plots <- list()
# maxvar <- family_variance_within_motus %>%
#     mutate(Genus = str_split_fixed(Genus, " ", n = 3)[, 2]) %>%
#     inner_join(data.frame(Genus = map_chr(names(mucin_pathway_colors), \(x) str_split(x, " ")[[1]][1]))) %>%
#     pull(var) %>%
#     max()
for (ge in 
    #c("Bacteroides")
    map_chr(names(mucin_pathway_colors), \(x) str_split(x, " ")[[1]][1])
) {
    print(str_c("Processing ", ge, "..."))
    now <- family_variance_within_motus %>%
        #filter(cazy_family != "GH2") %>%
        group_by(cazy_family) %>%
        filter(str_detect(Genus, ge))
    now <- now %>%
        group_by(cazy_family) %>%
        #filter(any(prevalence > 0.2)) %>% 
        filter(mean(prevalence > 0.5) > 0.2) %>% 
        select(mOTU_ID, Species, cazy_family, mean_copy_num, var, cov) %>%
        pivot_longer(c(var, mean_copy_num), names_to = "metric", values_to = "value") %>%
        identity() %>%
        mutate(metric = factor(metric, levels = c("mean_copy_num", "var"))) %>%
        mutate(frac = 0.5) %>%
        #inner_join(cazyAnnots %>% filter(GAG == "Yes" | Mucin == "Yes"), by = c("cazy_family" = "Subfamily")) %>%
        {
            if (only_mucin) {
                (.) %>% 
                    inner_join(data.frame(cazy_family = a)) %>%
                    mutate(cazy_family = factor(cazy_family, levels = a))
            } else {
                (.)
            }
        } %>%
        mutate(Species = str_replace(Species, "[0-9]+ ", "")) %>%
        mutate(Species = str_replace(Species, "NA ", "")) %>%
        mutate(Species = map_chr(Species, \(x) {
            xx <- str_split(x, " ")[[1]][1:2]
            return(str_c(xx, collapse = " "))
        })) %>%
        mutate(Species = str_replace(Species, "species", "sp."))

    unit_of_interest <- "mOTU_ID"
    column_name <- "Species"
    base_names <- now[[unit_of_interest]]
    column_info_to_add <- now[[column_name]]
    column_info_to_add <- str_replace_all(column_info_to_add, "\\[", "")
    column_info_to_add <- str_replace_all(column_info_to_add, "\\]", "")
    first_words <- map_chr(column_info_to_add, \(x) str_split(x, " ")[[1]][1])
    first_letters_with_dot <- str_c(str_sub(str_to_title(first_words), 1, 1), ". ")
    column_info_to_add <- str_replace(column_info_to_add, "[a-zA-Z]+ ", first_letters_with_dot)
    base_names <- str_c(column_info_to_add, " [", base_names, "]")
    now[[unit_of_interest]] <- base_names

    tmp <- now %>% 
        filter(metric == 'mean_copy_num') %>%
        pivot_wider(id_cols = mOTU_ID, names_from = cazy_family, values_from = value) %>% 
        #select(-mOTU_ID) %>% 
        column_to_rownames("mOTU_ID") %>%
        as.matrix()
    now <- now %>% 
        mutate(cazy_family = factor(cazy_family, levels = rev(colnames(tmp)[order(apply(tmp, 2, \(x) sum(x)))]))) %>%
        mutate(mOTU_ID = factor(mOTU_ID, levels = rownames(tmp)[order(apply(tmp, 1, \(x) sum(x)))]))
    w_tile <- 0.9
    h_tile <- 0.9

    heatmap_plot <- ggplot(now) +
    geom_split_tile(data = now[now$metric=="mean_copy_num", ] %>% mutate(mean_copy_num = value), aes(x = cazy_family, y = mOTU_ID, fill = log10(mean_copy_num+1), split = fct_rev(metric), frac=frac),colour = "white", linewidth = 0,width = w_tile, height = h_tile) +
    scale_fill_gradient(low = "white",high = "#1F77B4",na.value = "lightgrey", limits = c(0, 1.6)) +
    ggnewscale::new_scale_fill()+
    geom_split_tile(data = now[now$metric=="var", ] %>% mutate(variance = value) %>% mutate(variance = ifelse(variance > 15, 15, variance)), aes(x = cazy_family, y = mOTU_ID, fill = variance, split = fct_rev(metric), frac=frac),colour = "white", linewidth = 0,width = w_tile, height = h_tile) +
    scale_fill_gradient(low = "white",high = "#FF7F0E",na.value = "lightgrey", limits = c(0, 15)) +
    #scale_split(guide = guide_legend(override.aes = list(fill=c("#FFDAB9","#ADD8E6")))) +
    theme_presentation() +
    coord_fixed() +
    theme(
        #axis.text.x = element_text(angle = 45, hjust = 1,size = 8,face = ""),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
        )+
    scale_x_discrete(limits = a) +
    #ggtitle("MSI Tumors - HMR-seq vs TCGA WXS")+
    #labs(caption = "All associations with q>0.05 of MMR associated hallmark pathways and tumor enriched genera are shown.")+
    NULL
    # top_annot <- cazyAnnots %>%
    #     filter(GAG == "Yes" | Mucin == "Yes") %>%
    #     select(Subfamily, GAG, Mucin) %>%
    #     as_tibble()  %>%
    #     inner_join(data.frame(Subfamily = levels(now$cazy_family))) %>%
    #     mutate(Subfamily = factor(Subfamily, levels = hclust_o_cols$labels[hclust_o_cols$order])) %>%
    #     rename(cazy_family = Subfamily) %>%
    #     pivot_longer(-cazy_family, names_to = "annotation", values_to = "value") %>%
    #     ggplot(aes(x = cazy_family, y = annotation, value = value)) +
    #     theme_presentation() + 
    #     geom_tile(aes(fill = value), color = "white") +
    #     scale_fill_manual(values = c("white", "black")) + 
    #     theme(
    #         axis.text.x = element_blank(),
    #         axis.ticks.x = element_blank(),
    #         axis.title.x = element_blank(),
    #         axis.title.y = element_blank()) +
    library(patchwork)
    ggsave(plot = heatmap_plot, filename = here('figures', "revisions", str_c(ge, "_split_heatmap.pdf")), width = 8, height = 4)
    all_plots[[ge]] <- heatmap_plot
}

all_plots[1:(length(all_plots) - 1 )] <- map(all_plots[1:(length(all_plots) - 1 )], \(x) {
    return(
        x +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank()
            )
    )
})
ggsave(plot = wrap_plots(all_plots, ncol = 1) + plot_layout(guides = 'collect'), filename = here('figures', "revisions", "mucin_split_heatmap_abundances.pdf"), width = 8, height = 7.5)


now <- family_variance_within_motus %>%
    group_by(mOTU_ID, cazy_family) %>%
    filter(any(prevalence > 0.2)) %>% 
    group_by(cazy_family) %>%
    select(mOTU_ID, Genus, Species, cazy_family, var, cov, prevalence) %>%
    mutate(variance = var) %>%
    inner_join(cazyAnnots %>% filter(GAG == "Yes" | Mucin == "Yes"), by = c("cazy_family" = "Subfamily")) %>%
    mutate(Genus = str_replace(Genus, "[0-9]+ ", "")) %>%
    mutate(Species = str_replace(Species, "[0-9]+ ", "")) %>%
    mutate(Species = str_replace(Species, "NA ", "")) %>%
    mutate(Species = map_chr(Species, \(x) {
        xx <- str_split(x, " ")[[1]][1:2]
        return(str_c(xx, collapse = " "))
    })) %>%
    inner_join(data.frame(Genus = c(
        "Akkermansia",
        "Bacteroides",
        "Barnesiella",
        "Coprobacter",
        "Eisenbergiella",
        "Hungatella",
        "Parabacteroides",
        "Paraprevotella"
    ))) 
unit_of_interest <- "mOTU_ID"
column_name <- "Species"
base_names <- now[[unit_of_interest]]
column_info_to_add <- now[[column_name]]
column_info_to_add <- str_replace_all(column_info_to_add, "\\[", "")
column_info_to_add <- str_replace_all(column_info_to_add, "\\]", "")
first_words <- map_chr(column_info_to_add, \(x) str_split(x, " ")[[1]][1])
first_letters_with_dot <- str_c(str_sub(str_to_title(first_words), 1, 1), ". ")
column_info_to_add <- str_replace(column_info_to_add, "[a-zA-Z]+ ", first_letters_with_dot)
base_names <- str_c(column_info_to_add, " [", base_names, "]")
now[[unit_of_interest]] <- base_names

library(ggrepel)
dat <- now %>%
    filter(prevalence > 0.2) %>%
    inner_join(data.frame(cazy_family = a)) %>%
    mutate(cazy_family = factor(cazy_family, levels = a)) %>%
    select(mOTU_ID, Genus, Species, cazy_family, var, prevalence) %>%
    pivot_wider(names_from = c(cazy_family), values_from = c(var, prevalence), values_fill = 0) %>%
    pivot_longer(-c(mOTU_ID, Genus, Species), names_to = "cazy_family", values_to = "value") %>%
    mutate(metric = factor(str_split_fixed(cazy_family, "_", n = 3)[, 1], levels = c("var", "prevalence"))) %>%
    mutate(cazy_family = str_split_fixed(cazy_family, "_", n = 3)[, 2]) %>%
    pivot_wider(names_from = metric, values_from = c(value)) %>%
    mutate(cazy_family = factor(cazy_family, levels = a))
names(mucin_pathway_colors) <- map_chr(names(mucin_pathway_colors), \(x) str_split(x, " ")[[1]][1])
po_jd <- position_jitterdodge(jitter.width = 0.05, dodge.width = 0.75)
ggsave(
    plot = 
    ggplot(data = dat, aes(x = cazy_family, y = var, color = Genus, group = Genus)) +
    #geom_boxplot() +
    geom_point(position = po_jd, alpha = 0.5) +
    geom_text_repel(
        data = dat %>%
         filter(var > 0.15) %>%
         mutate(sp = map_chr(mOTU_ID, \(x) {
            tmp <- str_split(x, " ")[[1]][1:2]
            return(str_c(tmp, collapse = " "))
         })), aes(label = sp, color = Genus), size = 3, position = po_jd) +
    theme_presentation() +
    ylab("Variance (presence/absence)") +
    scale_color_manual(values = mucin_pathway_colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)), filename = here('figures', "revisions", "variance_beeswarm.pdf"), width = 12, height = 4)



    

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
# get_whole_deal_motu_level(
#     genus_number_and_name = "1263 Ruminococcus",
#     taxa_to_add_explicitly = c("torques", "gnavus", "obeum"), data_transformation = "zscore"
# )
# get_whole_deal_motu_level(
#     genus_number_and_name = "1263 Ruminococcus",
#     taxa_to_add_explicitly = c("torques", "gnavus", "obeum"), 
#     transform_to_prevalence = FALSE,
#     data_transformation = "log10"
# )
get_whole_deal_motu_level(
    genus_number_and_name = "1263 Ruminococcus",
    #taxa_to_add_explicitly = c("torques", "gnavus", "obeum"), 
    taxa_to_add_explicitly = NULL, 
    transform_to_prevalence = TRUE,
    data_transformation = "NONE"
)
# get_whole_deal_genome_level(
#     taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "gnavus")) %>% pull(mOTU_ID),
#     genus = "Ruminococcus",
#     data_transformation = 'zscore',
#     families_of_interest = NULL)
# get_whole_deal_genome_level(
#     taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "gnavus")) %>% pull(mOTU_ID),
#     genus = "Ruminococcus",
#     data_transformation = 'log10',
#     transform_to_prevalence = FALSE,
#     families_of_interest = NULL)
get_whole_deal_genome_level(
    taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "gnavus")) %>% pull(mOTU_ID),
    genus = "Ruminococcus",
    data_transformation = 'NONE',
    transform_to_prevalence = TRUE,
    families_of_interest = NULL)
# get_whole_deal_genome_level(
#     taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "torques")) %>% pull(mOTU_ID),
#     genus = "Ruminococcus",
#     data_transformation = 'zscore',
#     families_of_interest = NULL)
# get_whole_deal_genome_level(
#     taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "torques")) %>% pull(mOTU_ID),
#     genus = "Ruminococcus",
#     data_transformation = 'log10',
#     families_of_interest = NULL)
## We only have one species of Eisenbergiella, so we can't really do much here...
# get_whole_deal_motu_level(
#     genus_number_and_name = "1432051 Eisenbergiella",
#     taxa_to_add_explicitly = NULL
# )
source(here("scripts", "revisions_aditions", "explore_within_genus_cazy_diversity_utils.r"))
# get_whole_deal_genome_level(
#     taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "Eisenbergiella")) %>% pull(mOTU_ID),
#     genus = "Eisenbergiella",
#     data_transformation = 'zscore',
#     families_of_interest = NULL)
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
# get_whole_deal_genome_level(
#     taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "Hungatella")) %>% pull(mOTU_ID),
#     genus = "Hungatella",
#     data_transformation = 'zscore',
#     families_of_interest = NULL)
get_whole_deal_genome_level(
    taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "Hungatella")) %>% pull(mOTU_ID),
    genus = "Hungatella",
    data_transformation = 'log10',
    families_of_interest = NULL)
# get_whole_deal_motu_level(
#     genus_number_and_name = "816 Bacteroides",
#     data_transformation = 'zscore'
# )
get_whole_deal_motu_level(
    genus_number_and_name = "816 Bacteroides",
    data_transformation = 'log10'
)
# This runs for too long - heatmap is very large
# get_whole_deal_genome_level(
#     taxon_of_interest = motus_level_agg %>% filter(str_detect(Species, "Bacteroides")) %>% filter(number_genomes > 500) %>% pull(mOTU_ID),
#     genus = "Bacteroides",
#     families_of_interest = NULL)

set.seed(42)
almeida_paths <- read_tsv('/g/scb2/zeller/SHARED/DATA/assembled_genomes/Almeida_2020_combined_set/find.prokka.whole.path', col_names = F) %>% 
    mutate(genome = str_split_fixed(X1, "/", n= 15)[, 13]) %>% 
    mutate(genome = str_replace(genome, ".faa", "")) %>%
    rename(genome = genome, faa_path = X1)
# sample_conversion <- read_tsv('/g/scb/bork/fullam/proGenomes/freeze13/sample_conversion.tsv', col_names = F)
# pg3_represenatives <- read_tsv('/g/scb/bork/fullam/proGenomes/freeze13/specI_clusters/final_cluster_names/representatives.tsv', col_names = F) %>%
#     rename(pg3_gca = X2, genome = X3) %>%
#     select(pg3_gca, genome)
# pg3_paths <- read_tsv('/g/scb/bork/fullam/proGenomes/freeze13/all_files_metadata.tsv_v0') %>%
#     rename(pg3_gca = sample_id, pg3_path = genbank_path) %>%
#     select(pg3_gca, pg3_path)
representative_genomes_by_motu <- genome_level %>% 
    filter(mOTU_ID != "") %>% 
    group_by(mOTU_ID) %>% 
    filter(n() >= 50) %>%
    nest() %>%
    mutate(contains_isolate_genome = map_lgl(data, \(x) {
        if (any(x$Genome_type == "Isolate")) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    })) %>%
    mutate(almeida_genome = map2(data, contains_isolate_genome, \(x, cont) {
        if (cont) {
            return(x %>%
                filter(Genome_type == "Isolate") %>%
                sample_n(1) %>%
                select(genome, Genome_type))
        } else {
            me <- median(x$Length)
            return(x %>% 
                mutate(diff = abs(Length - me)) %>%
                arrange(diff) %>%
                head(1) %>%
                select(genome, Genome_type))
        }

    })) %>%
    select(-data) %>%
    unnest()

representative_genomes_by_motu %>%
    left_join(almeida_paths) %>%
    ungroup() %>%
    select(-c(mOTU_ID, Genome_type, )) %>%
    inner_join(genome_level, by = 'genome') %>%
    write_tsv(here('data', 'paths_to_276_rep_genomes.txt'))

##############################################################
# Check this here for the preliminary  final list of genomes 
##############################################################