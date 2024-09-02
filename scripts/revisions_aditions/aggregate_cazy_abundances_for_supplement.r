library(tidyverse)
library(here)
test_run <- FALSE
reload_data <- TRUE
if (reload_data) {
    cazy_all <- read_tsv('/g/scb/zeller/karcher/cayman_paper/data/almeida_cazy_annotations.tsv')
}

# Aggregate mean[sd] of CAZy genes per mOTUs level cluster (together with it's taxonomy)
motus_level_agg <- cazy_all %>%
    {
        if (test_run) {
            (.) %>% head(5000)
        } else {
            (.)
        }
    } %>%
    group_by(
        genome,
        mOTU_ID,
        Kingdom,
        Phylum,
        Class,
        Order,
        Family,
        Genus,
        Species,
        cazy_family) %>%
    summarize(
        copy_number = n()
        # copy_number = length(unique(sequenceID))
    ) %>%
    group_by(
        mOTU_ID,
        Kingdom,
        Phylum,
        Class,
        Order,
        Family,
        Genus,
        Species,
        cazy_family) %>%
    summarize(
        mean_copy_number = mean(copy_number),
        sd_copy_number = sd(copy_number)
    )

motus_level_agg %>%
    mutate(
        mean_copy_number = as.character(mean_copy_number),
        sd_copy_number = ifelse(is.na(sd_copy_number), "NA", as.character(sd_copy_number))) %>%
    mutate(cell_entries = str_c(mean_copy_number, ' [', sd_copy_number, ']')) %>%
    select(-mean_copy_number, -sd_copy_number) %>%
    pivot_wider(
        names_from = cazy_family,
        values_from = cell_entries,
        values_fill = list(cell_entries = '0 [NA]')) %>%
    filter(!is.na(mOTU_ID)) %>%
    left_join(cazy_all %>% select(genome, mOTU_ID) %>% ungroup() %>% distinct() %>% group_by(mOTU_ID) %>% tally() %>% select(mOTU_ID, n) %>% rename(number_genomes = n), by = 'mOTU_ID') %>%
    relocate(
        "mOTU_ID",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "number_genomes") %>%
    arrange(desc(number_genomes)) %>%
    write_tsv(here('data', 'motus_level_mean_cazy_abundances.tsv'))

## Same as above but species-level

# Aggregate mean[sd] of CAZy genes per mOTUs level cluster (together with it's taxonomy)
genome_level <- cazy_all %>%
    {
        if (test_run) {
            (.) %>% head(500000)
        } else {
            (.)
        }
    } %>%
    group_by(
        genome,
        mOTU_ID,
        Kingdom,
        Phylum,
        Class,
        Order,
        Family,
        Genus,
        Species,
        cazy_family) %>%
    tally() %>%
    rename(copy_number = n)

genome_level %>%
    pivot_wider(
        names_from = cazy_family,
        values_from = copy_number,
        values_fill = 0) %>%
    write_tsv(here('data', 'genome_level_cazy_abundances.tsv'))
