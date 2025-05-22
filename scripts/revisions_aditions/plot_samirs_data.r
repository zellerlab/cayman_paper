# conda activate /g/scb2/zeller/karcher/mambaforge/envs/r_growth_analysis
library(readxl)
library(tidyverse)
#library(pracma)
library(zoo)
#library(ggembl)

rna_selected_strains <- read_csv('/g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/selected_strains_for_RNA_seq.csv')

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
#data <- read_xlsx('/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run2.xlsx') %>%
#data <- read_xlsx('/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run3.xlsx') %>%
#data <- read_xlsx('/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run4_1-1.xlsx') %>%
# Based on email from Samir from February, Run4-2 should be correct for Eisenbergiella and Hungatella
data <- read_xlsx('/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run4_2-1.xlsx') %>%
    mutate(time_h = as.numeric(time)) %>%
    rename(well = Variable) %>%
    inner_join(data.frame(species = c("Eisenbergiella tayi", "Hungatella hathewayi"))) %>%
    # Remove 2 noisy measurements
    anti_join(
        tibble(
            plate = c("P1", "P1"),
            media = c("WCA", "WCA"),
            time = c(21, 23),
            well = c("E3", "F3")
        )
    ) %>%
    mutate(strain = str_c(species, "_", strainID)) %>%
    mutate(
        strain = case_when(
            strain == "Eisenbergiella tayi_DSM26961" ~ "Eisenbergiella tayi (DSM26961)",
            strain == "Hungatella hathewayi_DSM13479" ~ "Hungatella hathewayi (DSM13479)",
        )
    ) %>%
    mutate(
        condition = case_when(
            condition == "with_mucin_III_0.5%" ~ "With mucin",
            condition == "without_mucin" ~ "Without mucin"
        )
    ) %>%
    filter(time <= 30)

# plot mgam growth curves
for (medium in unique(data$media)) {
    p <- data %>%
        filter(media == medium) %>%
        mutate(g = str_c(well, plate)) %>%
        ggplot() +
        scale_fill_identity() +
        geom_line(aes(x = time_h, y = OD_adjusted, color = condition, group = g)) +
        theme_classic() +
        facet_wrap(. ~ strain, ncol = 3) +
        # make facet_wrap text size smaller
        theme(
            strip.text = element_text(size = 8)
        ) +
        xlab("time [h]") +
        ylab("OD") +
        scale_y_log10() +
        scale_color_manual(
            values = c(
                "With mucin" = "#8ecc52",
                "Without mucin" = "#808080"
            )
        ) +
        geom_vline(
            data = data.frame(
                xintercept = c(19, 20),
                strain = c("Eisenbergiella tayi (DSM26961)", "Hungatella hathewayi (DSM13479)")
            ),
            aes(xintercept = xintercept),
            alpha = 0.5,
            linetype = "longdash",
        )
        
        NULL

    #ggsave(plot = p, filename = "/g/scb/zeller/karcher/cayman_paper/figures/revisions/mgam_growth_curves.pdf", width = 9, height = 6)
    ggsave(plot = p, filename = str_c("/g/scb/zeller/karcher/cayman_paper/figures/revisions/", medium, "_growth_curves.pdf"), width = 6, height = 2.25)
}
