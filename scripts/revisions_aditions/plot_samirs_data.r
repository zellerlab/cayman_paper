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
        "Eisenbergiella tayi"
        ),
    purpose = c(
        "assayed", 
        "assayed", 
        "assayed", 
        "assayed", 
        "assayed", 
        "assayed")
) %>%
    distinct() %>%
    mutate(col = case_when(
        purpose == "assayed" ~ "#006400",
        purpose == "pos ctrl" ~ "#00008B",
        purpose == "neg ctrl" ~ "#808080",
        purpose == 'unclear' ~ "#FFFFFF"
    ))
df$col <- paste0(df$col, "40")

#data <- read_xlsx('/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run4_2-1.xlsx') %>%
data <- read_tsv('/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run4_2-1_subset.tsv') %>%
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
            strain == "Eisenbergiella tayi_DSM26961" ~ "E. tayi (DSM26961)",
            strain == "Hungatella hathewayi_DSM13479" ~ "H. hathewayi (DSM13479)",
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
        geom_vline(
            data = data.frame(
                xintercept = c(19, 20),
                strain = c("E. tayi (DSM26961)", "H. hathewayi (DSM13479)")
            ),
            aes(xintercept = xintercept),
            alpha = 0.5,
            linetype = "longdash",
        ) +        
        scale_fill_identity() +
        geom_line(aes(x = time_h, y = OD_adjusted, color = condition, group = g)) +
        theme_classic() +
        facet_wrap(. ~ strain, ncol = 1) +
        # make facet_wrap text size smaller
        theme(
            strip.text = element_text(size = 8)
        ) +
        xlab("time [h]") +
        ylab("OD") +
        scale_y_log10(
            breaks = c(0.01, 0.1, 1),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
            ) +
        scale_color_manual(
            values = c(
                "With mucin" = "#8ecc52",
                "Without mucin" = "#808080"
            )
        ) +
        # put legend to bottom in horizontal fasghion
        theme(
            axis.text.x = element_text(size = 7),
            axis.text.y = element_text(size = 7)
            #legend.position = "bottom",          # Place legend at the bottom
            #legend.direction = "horizontal",     # Make legend horizontal
            #legend.box = "vertical"      
        ) +     
        guides(color = guide_legend(ncol = 1)) +
        annotation_logticks(
            sides = 'l', 
            size = 0.25,
            short = unit(.5,"mm"),
            mid = unit(1.25,"mm"),
            long = unit(1.75,"mm")            
            )   +
        NULL
    #ggsave(plot = p, filename = "/g/scb/zeller/karcher/cayman_paper/figures/revisions/mgam_growth_curves.pdf", width = 9, height = 6)
    ggsave(plot = p, filename = str_c("/g/scb/zeller/karcher/cayman_paper/figures/revisions/", medium, "_growth_curves.pdf"), width = 3.8*0.9, height = 3.8*0.8)
}
