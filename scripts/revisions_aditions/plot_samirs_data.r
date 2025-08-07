# conda activate /g/scb2/zeller/karcher/mambaforge/envs/r_growth_analysis
library(readxl)
library(tidyverse)
#library(pracma)
library(zoo)
#library(ggembl)

rna_selected_strains <- read_csv('/g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/selected_strains_for_RNA_seq.csv')
startingOD <- 0.7 * 10^-4

interpolate_points <- function(data, n = 100) {
 
    # Ensure the data is sorted by x
    data <- data %>% arrange(time_h)
    
    # Perform interpolation
    interpolated <- approx(
        x = data$time_h,
        y = log10(data$OD + startingOD),
        n = n,  # Number of points to interpolate
        method = "linear"  # Linear interpolation
    )

    # Create a new data frame with the interpolated values for time_h and OD and all the other values in data as before
    clns <- colnames(data)
    clns <- clns[!clns %in% c("time_h", "OD")]

    # Replicate the other columns for each interpolated point
    other_columns <- data[1, clns, drop = FALSE]  # Take the first row's other columns
    other_columns <- other_columns[rep(1, length(interpolated$x)), ]  # Replicate for each interpolated point

    # Combine interpolated values with replicated columns
    interpolated_data <- cbind(
    time_h = interpolated$x,
    OD = interpolated$y,
    other_columns
    )
    return(interpolated_data)
}

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
    select(-OD_adjusted) %>%
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

da <- data_parsed_bound %>%
    anti_join(data.frame(media = "LBA")) %>%
    # Choose batches with smallest intra-batch variability over replicates
    inner_join(
        data.frame(
            media = c(
                "GMM+LAB",
                "SB",
                "WCA",
                "mGAM"
            ),
            batch = c(
                "1",
                "3",
                "3",
                "3"
            )
        )
    ) %>%
    filter(!(media == "WCA" & batch == "3" & strainID %in% c("DSM26961", "DSM13479"))) %>%
    rbind(
        data_parsed_bound %>%
            filter(media == "WCA" & batch == "5" & strainID %in% c("DSM26961", "DSM13479")) %>%
            # Remove faulty measurements
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

set.seed(1)
da <- da %>%
    group_by(species, media, condition, g) %>% nest() %>% group_by(species, media, condition) %>%
    sample_n(3) %>%
    unnest() %>%
    group_by(
        media
    )

da$medium_batch <- factor(da$medium_batch, levels = c('mGAM', sort(unique(da$medium_batch))[!sort(unique(da$medium_batch)) %in% c('mGAM')]))

daa <- da %>% 
    select(
        -c(time, plate, well, batch, replicate)
    ) %>%
    filter(OD > 0.01 | media != 'WCA') %>% 
    rbind(
        da %>% 
            filter(media == "WCA") %>% 
            select(strain, strainID, species, medium_batch, condition, g) %>% 
            distinct() %>% 
            mutate(OD = 0) %>% 
            mutate(time_h = 0)
    ) %>%
    mutate(growth_inferred = OD < 0.01) %>%
    filter(
        !(strain == "H. hathewayi\n(DSM13479)" & OD > 0.01 & OD < 0.015)
    )

get_x_at_y <- function(model, y_value) {
  # Extract coefficients from the linear model
  intercept <- coef(model)[1]  # Intercept (b)
  slope <- coef(model)[2]      # Slope (m)
  
  # Calculate x using the formula: x = (y - b) / m
  x_value <- (y_value - intercept) / slope
  
  return(x_value)
}

fitwindowsize <- 3
te <- daa %>%
    filter(media == "WCA") %>% filter(condition == "Without mucin") %>%
    filter(!growth_inferred) %>%
    group_by(condition, media, strainID,  medium_batch, species, strain, g)  %>%
    nest() %>%
    mutate(
        windows = map(data, \(x) {
            windowss <- list()
            for (startindex in c(1:10)) {
                tmp <- x[startindex:dim(x)[1], ] %>%
                    head(fitwindowsize)
                windowss[[length(windowss) + 1]] <- tmp
                names(windowss)[length(windowss)] <- str_c("window_", startindex)
            }
            return(windowss)
        }
        )
    ) %>%
    mutate(linfit = map(
        windows, \(x) {
            rsqs <- list()
            fits <- list()
            slopes <- list()
            for (window in x) {
                if (any(is.na(window))) {
                    fits[[length(fits) + 1]] <- NA
                    rsqs[[length(rsqs) + 1]] <- NA
                    slopes[[length(slopes) + 1]] <- NA
                    next
                }
                fit <- lm(log10(OD) ~ time_h, data = window)    
                #fits <- c(fits, fit)
                fits[[length(fits) + 1]] <- fit
                goodness_of_fit <- summary(fit)$r.squared
                if (is.na(goodness_of_fit)) {
                    goodness_of_fit <- 0
                }
                rsqs[[length(rsqs) + 1]] <- goodness_of_fit
                slopes[[length(slopes) + 1]] <- fit$coefficients[2]

            }
            if (length(rsqs) == 0) {
                return(NA)
            }
            if (!any(rsqs > 0.8)) {   
                return(NA)
            }
            tmp <- data.frame(
                i = 1:length(x), r = unlist(rsqs), s = unlist(slopes)
            ) %>%
                mutate(good = r>0.8) %>%
                filter(good) %>%
                arrange(desc(s))
            
            index_of_best_rqs <- tmp$i[1]
            #print(tmp)
            
            fit <- fits[[index_of_best_rqs]]
            print(fit)
            return(fit)
        }
    )) %>%
    filter(!is.na(linfit)) %>%
    mutate(
    slope = map_dbl(linfit, \(x) {
        return(x$coefficients[2])
    }),
    intercept = map_dbl(linfit, \(x) {
        return(x$coefficients[1])
    }),
    xintercept = map_dbl(linfit, \(x) {
        return(-x$coefficients[1] / x$coefficients[2])
    }),
    x_at_startOD = map_dbl(linfit, \(x) {
        return(get_x_at_y(x, log10(0.7*(10^-4))))
    })
    ) %>%
    identity() %>%
    select(condition, media,  strainID, species, strain, g, medium_batch, x_at_startOD, slope, intercept) %>%
    filter(slope > 0) %>%
    identity()


daa <- daa %>% 
    filter(!growth_inferred | media != 'WCA') %>% 
    rbind(
    daa %>% 
        filter(growth_inferred & media == "WCA") %>%
        ungroup() %>%
        select(strain, media, species, strainID, medium_batch, condition, g) %>%
        distinct() %>%
        mutate(OD = 0.7*(10^-4)) %>%
        mutate(time_h = 0) %>%
    group_by(strainID, species, strain, medium_batch, condition, g) %>%
    nest() %>%
    ungroup() %>%
    left_join(te, by = c('condition', 'strainID', "species", 'strain', "g", "medium_batch")) %>%
    mutate(data = map2(
        data, x_at_startOD, \(dat, x_at_start) {
            if (is.na(x_at_start)) {
                return(dat)
            }
            if (x_at_start < 0) {
                return(dat)
            } else {
                dat <- rbind(
                    dat, 
                    data.frame(
                        media = dat$media[1],
                        OD = 0.7*(10^-4),
                        time_h = x_at_start
                    )
                )
                return(dat)
            }
        }
    )) %>%
    unnest() %>%
    select(-x_at_startOD)
)

daa <- daa  %>%
    select(-growth_inferred) %>%
    group_by(
        media, condition, strainID, species, strain, g
    ) %>%
    nest() %>%
    mutate(data = map2(data,  media, \(x, me) {
        # if (me == "WCA") {
        #     return(x)
        # }
        x <- interpolate_points(x)
        return(x)
    })) %>%
    unnest() %>%
    # REverse the pseudocount and log-transform in interpolate_pooints
    mutate(OD = 10^OD - startingOD) %>%
    mutate(growth_inferred = OD < 0.01)

p <-  ggplot() +
        geom_abline(
            slope = 0, intercept = log10(0.01), linetype = 'dashed', alpha = 0.3
        ) +
        geom_line(data = da, aes(x = time_h, y = OD + 0.01, color = condition, group = g)) +
        theme_classic() +
        facet_grid(strain ~ medium_batch) +
        xlab("Time [h]") +
        ylab("OD") +
        scale_y_log10(
            breaks = c(0.001, 0.01, 0.1, 1),
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
        ) +     
        guides(color = guide_legend(ncol = 1)) +
        annotation_logticks(
            sides = 'l', 
            size = 0.25,
            short = unit(.5,"mm"),
            mid = unit(1.25,"mm"),
            long = unit(1.75,"mm")            
            )   +
        scale_fill_identity() +
        theme(
            strip.text.x = element_text(size = 7),
            strip.text.y = element_text(angle = 0, size = 7),
            plot.title = element_text(size = 20),  # Increase the font size (e.g., 20)
            legend.position = "bottom",          # Place legend at the bottom
            legend.direction = "horizontal"     # Make legend horizontal            
        ) +
        NULL

p_wca <-  ggplot() +
        geom_abline(
            slope = 0, intercept = log10(0.01), linetype = 'dashed', alpha = 0.3
        ) +
        geom_line(data = daa %>% filter(media == "WCA") %>% filter(growth_inferred), aes(x = time_h, y = OD + startingOD, color = condition, group = g), linetype = 'dashed', show.legend = FALSE) +
        geom_line(data = daa %>% filter(media == "WCA") %>% filter(!growth_inferred), aes(x = time_h, y = OD + startingOD, color = condition, group = g), show.legend = FALSE) +        
        theme_classic() +
        facet_grid(strain ~ medium_batch) +
        xlab("Time [h]") +
        ylab("OD") +
        scale_y_log10(
            breaks = c(0.0001, 0.0001, 0.001, 0.01, 0.1, 1),
            labels = scales::trans_format("log10", scales::math_format(10^.x)),
            limits = c(0.00005, 1.5)
            ) +
        scale_color_manual(
            values = c(
                "With mucin" = "#8ecc52",
                "Without mucin" = "#808080"
            )
        ) +
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
    plot = p + p_wca + plot_layout(widths = c(3, 0.79)) +
        plot_annotation(tag_levels = 'a'),  # This adds 'a' and 'b' to the panels
    filename = "/g/scb/zeller/karcher/cayman_paper/figures/revisions/condensed_growth_curves.pdf", 
    width = 8, 
    height = 4,
    units = "in"
)

flll <- data.frame(
    strain = c(
        "H. hathewayi\n(DSM13479)",
        "E. tayi\n(DSM26961)"
    ),
    xinterc = c(
        20, 
        18
    )
)

p <-  ggplot() +
        geom_abline(
            slope = 0, intercept = log10(0.01), linetype = 'dashed', alpha = 0.3
        ) +
        geom_line(data = daa %>% filter((str_detect(strain, "DSM26961") | str_detect(strain, "DSM13479")) & media == "WCA") %>% filter(growth_inferred), aes(x = time_h, y = OD + startingOD, color = condition, group = g), linetype = 'dashed') +
        geom_line(data = daa %>% filter((str_detect(strain, "DSM26961") | str_detect(strain, "DSM13479")) & media == "WCA") %>% filter(!growth_inferred), aes(x = time_h, y = OD + startingOD, color = condition, group = g)) +
        theme_classic() +
        facet_grid(strain ~ medium_batch) +
        xlab("Time [h]") +
        ylab("OD") +
        geom_vline(
            data = flll, 
            aes(xintercept = xinterc), 
            linetype = 'dashed', 
            color = 'black', 
            size = 0.3,
            alpha = 0.3
        ) +
        scale_y_log10(
            breaks = c(0.0001, 0.0001, 0.001, 0.01, 0.1, 1),
            labels = scales::trans_format("log10", scales::math_format(10^.x)),
            limits = c(0.00005, 1.5)
            ) +
        scale_color_manual(
            values = c(
                "With mucin" = "#8ecc52",
                "Without mucin" = "#808080"
            )
        ) +
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
    filename = "/g/scb/zeller/karcher/cayman_paper/figures/revisions/new_Fig3B.pdf", 
    width = 2.3*1.3, 
    height = 3*1.5,
    units = "in"
)

