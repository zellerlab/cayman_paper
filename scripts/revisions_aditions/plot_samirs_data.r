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

data_parsed <- list()
batchid <- 1
for (data_path in c(
	'/g/scb/zeller/karcher/cayman_paper/data/data_without_background_v2.xlsx',
	'/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run2.xlsx',
	'/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run3.xlsx',
	'/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run4_1-1.xlsx',
	'/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run4_2-1.xlsx'
)) {
	#data <- read_xlsx('/g/scb/zeller/karcher/cayman_paper/data/data_without_background_Run4_2-1.xlsx') %>%
	data <- read_xlsx(data_path) %>%
		mutate(time_h = as.numeric(time)) %>%
		rename(well = Variable) %>%
		mutate(
			species = map_chr(species, \(x) {
				first_word <- str_split(x, " ")[[1]][1]
				second_word <- str_split(x, " ")[[1]][2]
				first_letter_of_first_word <- str_split(first_word, "")[[1]][1]  %>% str_c(". ", sep = "")
				return(str_c(first_letter_of_first_word, second_word, collapse = ""))
			})
		) %>%
		mutate(strain = str_c(species, "\n(", strainID, ")")) %>%
		mutate(batch = batchid) %>%
		filter(time <= 30) %>%
		mutate(
			media = ifelse(media == "mGAM_0.5x", "mGAM", media),
			condition = ifelse(condition == "with_mucin_0.5%", "with_mucin_II_0.5%", condition)
		)

    # Fix the different OD values that were added by Samir
    data$OD <- data$OD_adjusted - min(data$OD_adjusted)

    # Add repliate information explicitly for first batch of experiments (since it's missing from file)
    if (batchid == 1) {
        data <- data %>% 
            group_by(media, plate, time, condition, strainID, species, strain) %>% 
            nest() %>% 
            mutate(data = map(data, \(x) x %>% mutate(replicate = 1:dim(.)[1]))) %>% 
            unnest()
    }

	data_parsed[[length(data_parsed) + 1 ]] <- data
 	batchid <- batchid + 1
}

data_parsed_bound <- bind_rows(data_parsed) %>%
    # We have already recomputed the original OD-values so we can now add a consistent pseudocount
    select(-OD_adjusted) %>%
    mutate(OD_adjusted = OD + 0.01) %>%    
    filter(
        str_detect(strain, "tayi") | 
        str_detect(strain, "hathew") |
        str_detect(strain, "mucini") |
        str_detect(strain, "secundus")
    ) %>%
	mutate(
		media = factor(media, levels = sort(unique(media))),
		condition = factor(condition, levels = sort(unique(condition))),
		#batch = factor(batch, levels = sort(unique(batch)))
        batch = factor(batch, levels = c(1,2,3,4,5))
	)
data_parsed_bound$strain <- factor(data_parsed_bound$strain, levels = sort(unique(data_parsed_bound$strain)))

pdf(file = "/g/scb/zeller/karcher/cayman_paper/figures/revisions/all_growth_curves.pdf", width = 12, height = 5)
for (medium in unique(data_parsed_bound$media)) {
    da <- data_parsed_bound %>%
        filter(media == medium) %>%
        mutate(g = str_c(well, plate)) %>%
        mutate(
            strain_sel = case_when(
                str_detect(strain, "DSM26961") ~ "#FFC0CB", 
                str_detect(strain, "13479") ~ "#FFC0CB", 
                TRUE ~ "#808080"
            )
        )
    p <- 
        ggplot(data = da) +
        # geom_vline(
        #     data = data.frame(
        #         xintercept = c(19, 20),
        #         strain = c("E. tayi (DSM26961)", "H. hathewayi (DSM13479)")
        #     ),
        #     aes(xintercept = xintercept),
        #     alpha = 0.5,
        #     linetype = "longdash",
        # ) +        
        geom_rect(
            data = da %>% ungroup() %>% select(species, strain, strain_sel) %>% distinct(),
            aes(fill = strain_sel),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            alpha = 0.3,
            show.legend = FALSE
        ) +
        scale_fill_identity() +
        geom_line(aes(x = time_h, y = OD_adjusted, color = condition, group = g)) +
        theme_classic() +
        facet_grid(strain ~ batch, drop = FALSE) +
        xlab("time [h]") +
        ylab("OD") +
        scale_y_log10(
            breaks = c(0.01, 0.1, 1),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
            ) +
        scale_color_manual(
            values = c(
                "with_mucin_III_0.5%" = "#8ecc52",
                "with_mucin_II_0.5%" = "#f2a900",
                "without_mucin" = "#808080"
            )
        ) +
        scale_fill_identity() + 
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
        ggtitle(medium) +
        theme(
            strip.text.x = element_text(size = 8),
            strip.text.y = element_text(angle = 0, size = 8),
            plot.title = element_text(size = 20)  # Increase the font size (e.g., 20)
        ) +
        NULL
    #ggsave(plot = p, filename = "/g/scb/zeller/karcher/cayman_paper/figures/revisions/mgam_growth_curves.pdf", width = 9, height = 6)
    ggsave(plot = p, filename = str_c("/g/scb/zeller/karcher/cayman_paper/figures/revisions/", medium, "_growth_curves.pdf"), width = 12, height = 7)
    print(p)
}
dev.off()

da <- data_parsed_bound %>%
    anti_join(data.frame(media = "LBA")) %>%
    inner_join(
        data.frame(
            media = c(
                "M17",
                "GMM+LAB",
                "SB",
                "WCA",
                "mGAM"
            ),
            batch = c(
                "2", 
                "1",
                "3",
                "3",
                "3"
            )
        )
    ) %>%
    # remove WCA growth data, batch 3, for DSM26961 and DSM13479... 
    filter(!(media == "WCA" & batch == "3" & strainID %in% c("DSM26961", "DSM13479"))) %>%
    # ... and replace with the data from batch 5
    rbind(
        data_parsed_bound %>%
            filter(media == "WCA" & batch == "5" & strainID %in% c("DSM26961", "DSM13479")) %>%
            # Remove faulty measuremnents
            anti_join(
                tibble(
                    plate = c("P1", "P1"),
                    media = c("WCA", "WCA"),
                    time = c(21, 23),
                    well = c("E3", "F3")
                )
            ) %>%
            mutate(
                batch = "3"
            )
    ) %>%
    inner_join(
        data.frame(strainID = c(
            # H. hathewayi
            "DSM13479",
            # A. muc
            "DSM22959",
            # E. tayi
            "DSM26961",
            # C. secundus
            "DSM28864/177"
        ))
    ) %>%
    mutate(g = str_c(well, plate)) %>%
    mutate(
        condition = factor(condition, levels = c(
            "without_mucin",
            "with_mucin_II_0.5%",
            "with_mucin_III_0.5%"
        )),
    ) %>%
    #filter(!condition == "with_mucin_III_0.5%") %>%
    mutate(condition = ifelse(
        str_detect(condition, "with_mucin"),
        "With mucin",
        "Without mucin"
    )) %>%
    group_by(
        media
    ) %>%
    nest() %>%
    mutate(data = map(data, \(x) {
        x$batch = as.numeric(as.factor(x$batch))
        return(x)
    })) %>%
    unnest() %>%
    mutate(medium_batch = media)

da$medium_batch <- factor(da$medium_batch, levels = c('mGAM', sort(unique(da$medium_batch))[!sort(unique(da$medium_batch)) %in% c('mGAM')]))

p <-   ggplot(data = da) +
        geom_line(aes(x = time_h, y = OD_adjusted, color = condition, group = g)) +
        theme_classic() +
        facet_grid(strain ~ medium_batch) +
        xlab("Time [h]") +
        ylab("OD") +
        scale_y_log10(
            breaks = c(0.01, 0.1, 1),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
            ) +
        scale_color_manual(
            values = c(
                #"with_mucin_III_0.5%" = "#8ecc52",
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
        #ggtitle(medium) +
        geom_rect(
            data = da %>%
                ungroup() %>%
                select(species, media, medium_batch, strain) %>%
                distinct() %>%
                mutate(col = ifelse(
                    (str_detect(strain, "DSM26961") | str_detect(strain, "DSM13479")) & media == "WCA",
                    "#FFC0CB40",
                    NA
                )),
                # Make it such that the background of each tile is colored
            aes(fill = col),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
        ) +
        scale_fill_identity() +
        theme(
            strip.text.x = element_text(size = 7),
            strip.text.y = element_text(angle = 0, size = 7),
            plot.title = element_text(size = 20),  # Increase the font size (e.g., 20)
            legend.position = "bottom",          # Place legend at the bottom
            legend.direction = "horizontal"     # Make legend horizontal            
            # put legend at bottom
        ) +
        NULL

ggsave(
    plot = p, 
    filename = "/g/scb/zeller/karcher/cayman_paper/figures/revisions/condensed_growth_curves.pdf", 
    width = 6, 
    height = 4,
    units = "in"
)
