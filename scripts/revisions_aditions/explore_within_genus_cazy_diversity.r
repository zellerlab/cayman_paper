##############
##############
### ATTENTION # Run from conda activate /g/scb/zeller/fspringe/Software/miniconda/envs/r_env_4.3.3
##############
##############

`library(tidyverse)
library(here)
library(readxl)
# library(pheatmap)
library(ComplexHeatmap)
library(ggembl)
library(ggnewscale)
library(patchwork)`
source(here('scripts', 'utils.r'))

pseudocount <- 1

# Mucin families shown in (old) Fig 2B, should now be Fig3A?
a <- families_2B <- c("GH2", "GH92", "GH20", "GH31", "GH29", "GH97", "GH95", "CBM32", "GH36", "GH35", 
  "GH33", "GH42", "GH130", "GH18", "GH16", "GH109", "CBM51", "GH110", "GH84", 
  "GH27", "GH89", "GH123", "GH85", "GH136", "GH112")

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
        # var = var(copy_number),
        # cov = sd(copy_number) / mean(copy_number)) %>%
        var = var(present),
        cov = sd(present) / mean(present)) %>%        
    filter(num_genomes > 10)

only_mucin <- TRUE
all_plots <- list()
for (ge in 
    #c("Bacteroides")
    map_chr(names(mucin_pathway_colors), \(x) str_split(x, " ")[[1]][1])
) {
    print(str_c("Processing ", ge, "..."))
    now <- family_variance_within_motus %>%
        #filter(cazy_family != "GH2") %>%
        group_by(cazy_family) %>%
        filter(str_detect(Genus, ge)) %>%
        group_by(cazy_family) %>%
        #filter(any(prevalence > 0.2)) %>% 
        filter(mean(prevalence > 0.5) > 0.2) %>% 
        select(mOTU_ID, Species, cazy_family, prevalence, var, cov) %>%
        pivot_longer(c(prevalence, var), names_to = "metric", values_to = "value") %>%
        identity() %>%
        mutate(metric = factor(metric, levels = c("prevalence", "var"))) %>%
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
        filter(metric == 'prevalence') %>%
        pivot_wider(id_cols = mOTU_ID, names_from = cazy_family, values_from = value) %>% 
        #select(-mOTU_ID) %>% 
        column_to_rownames("mOTU_ID") %>%
        as.matrix()
    now <- now %>% 
        mutate(cazy_family = factor(cazy_family, levels = rev(colnames(tmp)[order(apply(tmp, 2, \(x) sum(x)))]))) %>%
        mutate(mOTU_ID = factor(mOTU_ID, levels = rownames(tmp)[order(apply(tmp, 1, \(x) sum(x)))]))
    w_tile <- 0.9
    h_tile <- 0.9

    print(sum(is.na(now)))

    heatmap_plot <- ggplot(now) +
    geom_tile(data = now[now$metric=="prevalence", ] %>% mutate(prevalence = value), aes(x = cazy_family, y = mOTU_ID, fill = prevalence)) +
    theme_presentation() +
    coord_fixed(.65) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 7),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",    # Center the legend horizontally
        legend.box = "horizontal"                     
        )+
    scale_x_discrete(limits = a) +
    scale_fill_continuous(low = "white", high = "#1F77B4", na.value = "lightgrey", limits = c(0, 1)) +
    NULL

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
            axis.title.x = element_blank(),
            legend.position = "bottom",          # Place legend at the bottom
            legend.direction = "horizontal",    # Make legend horizontal
            legend.justification = "center"     # Center the legend            
            )
    )
})
ggsave(plot = wrap_plots(all_plots, ncol = 1) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
, filename = here('figures', "revisions", "mucin_split_heatmap.pdf"), width = 5, height = 7.5)

