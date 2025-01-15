# conda activate /g/scb2/zeller/karcher/mambaforge/envs/r_growth_analysis
library(readxl)
library(tidyverse)
#library(pracma)
library(zoo)
#library(ggembl)

get_auc_v1 <- function(x, y) {
    id <- order(x)
    AUC <- sum(diff(x[id])*rollmean(y[id],2))
    return(AUC)
}

get_auc_v2 <- function(x, y) { 
    AUC <- pracma::trapz(x, y)
    return(AUC)
}

df <- data.frame(
    species = c(
        "Hungatella hathewayi", 
        "Hungatella hathewayi", 
        "Hungatella hathewayi", 
        "Eisenbergiella tayi", 
        "Eisenbergiella tayi", 
        "Eisenbergiella tayi", 
        "Eisenbergiella massiliensis", 
        "Coprobacter fastidiosus", 
        "Coprobacter secundus", 
        "Akkermansia muciniphila", 
        "Pseudoflavonifractor capillosus", 
        "Eubacterium eligens", 
        "Roseburia intestinalis", 
        "Roseburia inulinivorans", 
        "Dorea formicigenerans", 
        "Dorea longicatena",
        "Clostridioides difficile",
        "Bacteroides uniformis",
        "Phocaeicola vulgatus"),
    purpose = c(
        "assayed", 
        "assayed", 
        "assayed", 
        "assayed", 
        "assayed", 
        "assayed", 
        "assayed", 
        "assayed", 
        "assayed", 
        "pos ctrl", 
        "neg ctrl", 
        "neg ctrl", 
        "neg ctrl", 
        "neg ctrl", 
        "neg ctrl", 
        "neg ctrl",
        "neg ctrl",
        'unclear',
        'unclear')
) %>%
    distinct() %>%
    mutate(col = case_when(
        purpose == "assayed" ~ "#006400",
        purpose == "pos ctrl" ~ "#00008B",
        purpose == "neg ctrl" ~ "#808080",
        purpose == 'unclear' ~ "#FFFFFF"
    ))
df$col <- paste0(df$col, "40")

#data <- read_xlsx('/g/scb/zeller/karcher/cayman_paper/data/data_without_background.xlsx') %>%
#data <- read_xlsx('/g/scb/zeller/karcher/cayman_paper/data/data_without_background_v2.xlsx') %>%
data <- read_xlsx('/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run2.xlsx') %>%
#data <- read_xlsx('/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run3.xlsx') %>%
    mutate(time_h = as.numeric(time)) %>%
    rename(well = Variable)

data %>%
    group_by(plate, well, media, time_h, condition, strainID) %>%
    tally()

auc_data <- data %>% 
    group_by(plate, well, media, condition, strainID, species) %>%
    summarize(
        auc = get_auc_v1(time_h, OD_adjusted),
        #auc = get_auc_v2(time_h, OD_adjusted)
    )

auc_data_agg <- auc_data %>%
    group_by(media, condition, strainID, species) %>%
    summarize(
        mean_auc = mean(auc), sd_auc = sd(auc)
        ) %>%
    left_join(df, by = 'species') %>%
    mutate(species = factor(species, levels = df$species))

species_values <- 1:length(unique(auc_data_agg$species))
pp <- position_dodge(width = 0.8)
p <- ggplot(auc_data_agg) +
    geom_tile(
        data = auc_data_agg %>%
            ungroup() %>% 
            select(species, col) %>%
            distinct(), aes(x = species, y = 0, fill = col), width = 1, height = Inf) +
    geom_point(aes(x = species, y = mean_auc, color = condition, group = strainID), position = pp) +
    geom_errorbar(aes(x = species, ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc, group = strainID), position = pp, width = 0.2, alpha = 0.5) +
    geom_vline(xintercept = species_values + 0.5, linetype = "dashed", color = "black") +
    scale_fill_identity() +
    theme_classic() +
    facet_wrap(media ~., ncol = 3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("AUC")
    NULL

ggsave(plot = p, filename = "/g/scb/zeller/karcher/cayman_paper/figures/revisions/growth_data_auc_v1.pdf", width = 10, height = 6)
#ggsave(plot = p, filename = "/g/scb/zeller/karcher/cayman_paper/figures/revisions/growth_data_auc_v2.pdf", width = 10, height = 6)

# plot mgam growth curves
for (medium in unique(data$media)) {
    p <- data %>%
        #filter(media == 'mGAM') %>%
        filter(media == medium) %>%
        mutate(strain = str_c(species, "_", strainID)) %>%
        ggplot() +
        geom_tile(
            data = auc_data_agg %>%
                ungroup() %>% 
                select(species, strainID, col) %>%
                mutate(strain = str_c(species, "_", strainID)) %>%
                distinct(), aes(x = 1, y = 0, fill = col), width = Inf, height = Inf) +    
        scale_fill_identity() +
        geom_line(aes(x = time_h, y = OD_adjusted, color = condition, group = well)) +
        theme_classic() +
        facet_wrap(. ~ strain, ncol = 3) +
        # make facet_wrap text size smaller
        theme(
            strip.text = element_text(size = 8)
        ) +
        xlab("time [h]") +
        scale_y_log10() +
        ylab("OD")
        NULL

    #ggsave(plot = p, filename = "/g/scb/zeller/karcher/cayman_paper/figures/revisions/mgam_growth_curves.pdf", width = 9, height = 6)
    ggsave(plot = p, filename = str_c("/g/scb/zeller/karcher/cayman_paper/figures/revisions/", medium, "_growth_curves.pdf"), width = 9, height = 6)
}

