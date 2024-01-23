library(here)
library(tidyverse)
library(ggsignif)
source(here('scripts', 'utils.r'))

p <- read_tsv(here('data', 'w_nw_profiles.tsv'))

(p %>%
    filter(genus == "g__Bacteroides") %>%
    group_by(sampleID, genus, non_westernized) %>%
    summarize(relAb = sum(relAb)) %>%
    rename(`Westernization status` = non_westernized) %>%
    mutate(`Westernization status` = ifelse(`Westernization status` == "Westernized", "Western", "Non-Western")) %>%
    ggplot(aes(x = `Westernization status`, y = relAb + 1E-5, fill = `Westernization status`)) +
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(trans = 'log10', breaks = c(1, 0.1, 0.01, 0.001, 0.0001, 0.000001)) + 
    scale_fill_manual(values = c("Western" = "#3B6FB6", 'Non-Western' = "#D41645")) +
    geom_signif(
    comparisons = list(c("Western", "Non-Western"))
    ) +
    ylab("Bacteroides\nrelative abundance"))  %>%
    ggsave(filename = here('figures', "Extended_Data_Fig6_Bacteroides_W_NW.pdf"), width = 5, height = 3.5)
