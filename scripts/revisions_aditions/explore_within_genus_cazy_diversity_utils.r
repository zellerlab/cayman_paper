save_heatmap <- function(o, path, dpi = 300, width = 8, height = 8) {
    pdf(file = path, width = width, height = height, pointsize = 12, family = "Helvetica")
    print(o)
    dev.off()
}

# breakpoint <- function() {
#     browser()
# }

get_base_col_name_boolean <- function(base_col_names) {
    return((!str_detect(base_col_names, "^CE") &
        !str_detect(base_col_names, "^CBM") &
        !str_detect(base_col_names, "^GH") &
        !str_detect(base_col_names, "^GT") &
        !str_detect(base_col_names, "^AA") &
        !str_detect(base_col_names, "^PL")))
}



get_whole_deal_motu_level <- function(
    genus_number_and_name,
    data_transformation,
    taxa_to_add_explicitly = NULL) {
    spli <- str_split(genus_number_and_name, " ")[[1]]
    stopifnot(length(spli) == 2)
    genus_name <- spli[2]
    genus_name_all_lowercase <- tolower(genus_name)
    genus_name_first_letter_uppercase <- str_to_title(genus_name_all_lowercase)
    dir.create(here("figures", "revisions", str_c(genus_name_first_letter_uppercase, '_cazy_zoom_in')), showWarnings = FALSE)
    if (!is.null(taxa_to_add_explicitly)) {
        taxa_to_add_explicitly <- str_c(taxa_to_add_explicitly, collapse = "|")
        taxa_to_add_explicitly <- motus_level_agg %>% filter(
            str_detect(Species, taxa_to_add_explicitly)) %>% pull(mOTU_ID) %>% map(\(x) c(x, "mOTU_ID"))
    }
    # pheatmap_object <- get_heatmap(
    #     data = motus_level_agg,
    #     taxon_of_interest = c(genus_number_and_name),
    #     level_of_interest = "Genus",
    #     unit_of_interest = 'mOTU_ID',
    #     data_transformation = 'zscore',
    #     families_of_interest = NULL,
    #     # taxa_to_add_explicitly = list(
    #     #     c("33038 [Ruminococcus] gnavus [[Ruminococcus] gnavus CAG:126/Clostridia bacterium UC5.1-2D9/[Ruminococcus] gnavus/Lachnospiraceae bacterium 2_1_58FAA]", "Species")),
    #     taxa_to_add_explicitly = taxa_to_add_explicitly,
    #     info_to_add = c("Species_simple"),
    #     column_annotation_names = c(
    #         "log10_#_genomes",
    #         "Family"
    #     )
    # )
    # save_heatmap(
    #     pheatmap_object,
    #     here("figures", "revisions", str_c(genus_name_first_letter_uppercase'_cazy_zoom_in'), str_c(genus_name_all_lowercase,"___all_families___heatmap.pdf")),
    #     width = 10, height = 4)
    pheatmap_object <- get_heatmap(
        data = motus_level_agg,
        taxon_of_interest = c(genus_number_and_name),
        level_of_interest = "Genus",
        unit_of_interest = 'mOTU_ID',
        data_transformation = data_transformation,
        families_of_interest = NULL,
        # taxa_to_add_explicitly = list(
        #     c("33038 [Ruminococcus] gnavus [[Ruminococcus] gnavus CAG:126/Clostridia bacterium UC5.1-2D9/[Ruminococcus] gnavus/Lachnospiraceae bacterium 2_1_58FAA]", "Species")),
        taxa_to_add_explicitly = taxa_to_add_explicitly,
        info_to_add = c("Species_simple"),
        column_annotation_names = c(
            "log10_#_genomes",
            # str_c("mean_dist_", data_transformation),
            # str_c("var_", data_transformation),
            "mean_dist",
            "var",
            "sum_copyn",
            "sum_uniq_f",
            # "sum_unique_cazyme_families",
            # "mean_dist_to_other_motus",
            # "variance_cazyme_z_scores",
            # "sum_cazyme_copy_number",
            # "sum_cazyme_families",
            ## "sum_unique_cazyme_families",
            "Family"
        )
    )
    save_heatmap(
        pheatmap_object,
        here("figures", "revisions", str_c(genus_name_first_letter_uppercase, '_cazy_zoom_in'), str_c(genus_name_all_lowercase, "data_transformation_", data_transformation, "___agnostic_family_selection___heatmap.pdf")),
        width = 16, height = 7)

    pheatmap_object <- get_heatmap(
        data = motus_level_agg,
        taxon_of_interest = c(genus_number_and_name),
        level_of_interest = "Genus",
        unit_of_interest = 'mOTU_ID',
        data_transformation = data_transformation,
        families_of_interest = cazyAnnots %>% filter(GAG == "Yes" | Mucin == "Yes") %>% pull(Subfamily) %>% unique(),
        # taxa_to_add_explicitly = list(
        #     c("33038 [Ruminococcus] gnavus [[Ruminococcus] gnavus CAG:126/Clostridia bacterium UC5.1-2D9/[Ruminococcus] gnavus/Lachnospiraceae bacterium 2_1_58FAA]", "Species")),
        taxa_to_add_explicitly = taxa_to_add_explicitly,
        info_to_add = c("Species_simple"),
        column_annotation_names = c(
            "log10_#_genomes",
            # str_c("mean_dist_", data_transformation),
            # str_c("var_", data_transformation),
            "mean_dist",
            "var",
            "sum_copyn",
            "sum_uniq_f",
            # "sum_unique_cazyme_families",
            # "mean_dist_to_other_motus",
            # "variance_cazyme_z_scores",
            # "sum_cazyme_copy_number",
            # "sum_cazyme_families",
            ## "sum_unique_cazyme_families",
            "Family"
        )
    )
    save_heatmap(
        pheatmap_object,
        here("figures", "revisions", str_c(genus_name_first_letter_uppercase, '_cazy_zoom_in'), str_c(genus_name_all_lowercase, "data_transformation_", data_transformation, "___mucin_gag___heatmap.pdf")),
        width = 16, height = 7)

    pheatmap_object <- get_heatmap(
        data = motus_level_agg,
        taxon_of_interest = c(genus_number_and_name),
        level_of_interest = "Genus",
        unit_of_interest = 'mOTU_ID',
        data_transformation = data_transformation,
        families_of_interest = c("GH95", "GH29"),
        # taxa_to_add_explicitly = list(
        #     c("33038 [Ruminococcus] gnavus [[Ruminococcus] gnavus CAG:126/Clostridia bacterium UC5.1-2D9/[Ruminococcus] gnavus/Lachnospiraceae bacterium 2_1_58FAA]", "Species")),
        taxa_to_add_explicitly = taxa_to_add_explicitly,
        info_to_add = c("Species_simple"),
        column_annotation_names = c(
            "log10_#_genomes",
            # str_c("mean_dist_", data_transformation),
            # str_c("var_", data_transformation),
            "mean_dist",
            "var",
            "sum_copyn",
            "sum_uniq_f",
            # "sum_unique_cazyme_families",
            # "mean_dist_to_other_motus",
            # "variance_cazyme_z_scores",
            # "sum_cazyme_copy_number",
            # "sum_cazyme_families",
            ## "sum_unique_cazyme_families",
            "Family"
        )
    )
    save_heatmap(
        pheatmap_object,
        here("figures", "revisions", str_c(genus_name_first_letter_uppercase, '_cazy_zoom_in'), str_c(genus_name_all_lowercase, "data_transformation_", data_transformation, "___GH29_GH95___heatmap.pdf")),
        width = 16, height = 7)

    pheatmap_object <- get_heatmap(
        data = motus_level_agg,
        taxon_of_interest = c(genus_number_and_name),
        level_of_interest = "Genus",
        unit_of_interest = 'mOTU_ID',
        data_transformation = data_transformation,
        families_of_interest = 'all_families',
        # taxa_to_add_explicitly = list(
        #     c("33038 [Ruminococcus] gnavus [[Ruminococcus] gnavus CAG:126/Clostridia bacterium UC5.1-2D9/[Ruminococcus] gnavus/Lachnospiraceae bacterium 2_1_58FAA]", "Species")),
        taxa_to_add_explicitly = taxa_to_add_explicitly,
        info_to_add = c("Species_simple"),
        column_annotation_names = c(
            "log10_#_genomes",
            # str_c("mean_dist_", data_transformation),
            # str_c("var_", data_transformation),
            "mean_dist",
            "var",
            "sum_copyn",
            "sum_uniq_f",
            # "sum_unique_cazyme_families",
            # "mean_dist_to_other_motus",
            # "variance_cazyme_z_scores",
            # "sum_cazyme_copy_number",
            # "sum_cazyme_families",
            ## "sum_unique_cazyme_families",
            "Family"
        )
    )
    save_heatmap(
        pheatmap_object,
        here("figures", "revisions", str_c(genus_name_first_letter_uppercase, '_cazy_zoom_in'), str_c(genus_name_all_lowercase, "data_transformation_", data_transformation, "___all_families___heatmap.pdf")),
        width = 16, height = 7)
}
get_whole_deal_genome_level <- function(
    taxon_of_interest,
    genus,
    data_transformation,
    families_of_interest = NULL) {
    pheatmap_object <- get_heatmap(
        data = genome_level,
        taxon_of_interest = taxon_of_interest,
        level_of_interest = "mOTU_ID",
        unit_of_interest = 'genome',
        data_transformation = data_transformation,
        families_of_interest = families_of_interest,
        column_annotation_names = c(
            "mOTU_ID"
        )
        # taxa_to_add_explicitly = list(c("33038 [Ruminococcus] gnavus [[Ruminococcus] gnavus CAG:126/Clostridia bacterium UC5.1-2D9/[Ruminococcus] gnavus/Lachnospiraceae bacterium 2_1_58FAA]", "Species")),
        # info_to_add = c("Species_simple")
        # info_to_add = NULL
    )
    save_heatmap(
        pheatmap_object,
        here("figures", "revisions", str_c(genus, '_cazy_zoom_in'), str_c(str_c(taxon_of_interest, collapse = "_"), "data_transformation_", data_transformation, "___agnostic_family_selection___heatmap.pdf")),
        width = 16, height = 7)

    pheatmap_object <- get_heatmap(
        data = genome_level,
        taxon_of_interest = taxon_of_interest,
        level_of_interest = "mOTU_ID",
        unit_of_interest = 'genome',
        data_transformation = data_transformation,
        families_of_interest = cazyAnnots %>% filter(GAG == "Yes" | Mucin == "Yes") %>% pull(Subfamily) %>% unique(),
        column_annotation_names = c(
            "mOTU_ID"
        )
        # taxa_to_add_explicitly = list(c("33038 [Ruminococcus] gnavus [[Ruminococcus] gnavus CAG:126/Clostridia bacterium UC5.1-2D9/[Ruminococcus] gnavus/Lachnospiraceae bacterium 2_1_58FAA]", "Species")),
        # info_to_add = c("Species_simple")
        # info_to_add = NULL
    )
    save_heatmap(
        pheatmap_object,
        here("figures", "revisions", str_c(genus, '_cazy_zoom_in'), str_c(str_c(taxon_of_interest, collapse = "_"), "data_transformation_", data_transformation, "___mucin_gag___heatmap.pdf")),
        width = 16, height = 7)

    pheatmap_object <- get_heatmap(
        data = genome_level,
        taxon_of_interest = taxon_of_interest,
        level_of_interest = "mOTU_ID",
        unit_of_interest = 'genome',
        data_transformation = data_transformation,
        families_of_interest = c("GH95", "GH29"),
        column_annotation_names = c(
            "mOTU_ID"
        )
        # taxa_to_add_explicitly = list(c("33038 [Ruminococcus] gnavus [[Ruminococcus] gnavus CAG:126/Clostridia bacterium UC5.1-2D9/[Ruminococcus] gnavus/Lachnospiraceae bacterium 2_1_58FAA]", "Species")),
        # info_to_add = c("Species_simple")
        # info_to_add = NULL
    )
    save_heatmap(
        pheatmap_object,
        here("figures", "revisions", str_c(genus, '_cazy_zoom_in'), str_c(str_c(taxon_of_interest, collapse = "_"), "data_transformation_", data_transformation, "___GH29_GH95___heatmap.pdf")),
        width = 16, height = 7)

    pheatmap_object <- get_heatmap(
        data = genome_level,
        taxon_of_interest = taxon_of_interest,
        level_of_interest = "mOTU_ID",
        unit_of_interest = 'genome',
        data_transformation = data_transformation,
        families_of_interest = 'all_families',
        column_annotation_names = c(
            "mOTU_ID"
        )
        # taxa_to_add_explicitly = list(c("33038 [Ruminococcus] gnavus [[Ruminococcus] gnavus CAG:126/Clostridia bacterium UC5.1-2D9/[Ruminococcus] gnavus/Lachnospiraceae bacterium 2_1_58FAA]", "Species")),
        # info_to_add = c("Species_simple")
        # info_to_add = NULL
    )
    save_heatmap(
        pheatmap_object,
        here("figures", "revisions", str_c(genus, '_cazy_zoom_in'), str_c(str_c(taxon_of_interest, collapse = "_"), "data_transformation_", data_transformation, "___all_families___heatmap.pdf")),
        width = 16, height = 7)
}

get_heatmap <- function(
    data,
    taxon_of_interest,
    level_of_interest,
    unit_of_interest,
    families_of_interest = NULL,
    taxa_to_add_explicitly = NULL,
    data_transformation = NULL,
    info_to_add = NULL,
    column_annotation_names = NULL
    ) {
    data_orig <- data
    print(str_c("Processing data for ", taxon_of_interest, " at ", level_of_interest, "..."))
    ij <- data.frame(tmp = taxon_of_interest)
    colnames(ij)[1] <- level_of_interest
    data <- data %>%
        inner_join(ij, by = level_of_interest)
    if (!is.null(taxa_to_add_explicitly)) {
        print("Adding taxa explicitly...")
        for (rowI in 1:length(taxa_to_add_explicitly)) {
            taxon <- taxa_to_add_explicitly[[rowI]][1]
            level <- taxa_to_add_explicitly[[rowI]][2]
            ij <- data.frame(tmp = taxon)
            colnames(ij)[1] <- level
            di <- dim(data)[1]
            data <- rbind(
                data,
                data_orig %>%
                    inner_join(ij))
            stopifnot(dim(data)[1] == (di + 1))
        }
    }
    data <- data %>% as_tibble()
    # if we're looking at aggregated data, only look at groups with at least 10 genomes
    # if ("number_genomes" %in% colnames(data)) {
    #     print("Removing groups with fewer than 20 genomes...")
    #     data <- data %>%
    #         filter(number_genomes >= 20)
    # }
    print(str_c("Retaining ", nrow(data), " groups..."))
    if (!is.null(families_of_interest)) {
        if (families_of_interest == "all_families") {
            print("Retaining all families...")
        } else {
            print(str_c("Filtering for ", length(families_of_interest), " families..."))
            base_col_names <- colnames(data)
            base_col_names <- base_col_names[
                get_base_col_name_boolean(base_col_names)]
            # data <- data[, c(base_col_names, families_of_interest)]
            data <- data[, colnames(data) %in% c(base_col_names, families_of_interest)]
            # Make sure all families are retained...
            if (dim(data)[2] == (length(base_col_names) + length(families_of_interest))) {
                print("All families retained.")
            } else {
                print("Some families were not retained...")
            }
        }

    } else {
        print(str_c("No families of interest provided, retaining all reasonably prevalent families..."))
        num_data <- data[, !get_base_col_name_boolean(colnames(data))]
        tmp <- apply(num_data, 2, \(x) mean(x >= 1) > 0.2)
        abundant_families <- names(tmp)[tmp]
        data <- data[, c(colnames(data)[get_base_col_name_boolean(colnames(data))], abundant_families)]
    }
    print("cleaning up dataframe...")
    if (!is.null(info_to_add)) {
        base_names <- data[[unit_of_interest]]
        for (column_name in info_to_add) {
            column_info_to_add <- data[[column_name]]
            column_info_to_add <- str_replace_all(column_info_to_add, "\\[", "")
            column_info_to_add <- str_replace_all(column_info_to_add, "\\]", "")
            first_words <- map_chr(column_info_to_add, \(x) str_split(x, " ")[[1]][1])
            first_letters_with_dot <- str_c(str_sub(str_to_title(first_words), 1, 1), ". ")
            column_info_to_add <- str_replace(column_info_to_add, "[a-zA-Z]+ ", first_letters_with_dot)
            base_names <- str_c(column_info_to_add, " [", base_names, "]")
        }
        data[[unit_of_interest]] <- base_names
    }

    data_all <- data
    data <- data %>%
        select(all_of(c(unit_of_interest, colnames(.)[!get_base_col_name_boolean(colnames(.))]))) %>%
        as.data.frame() %>%
        column_to_rownames(unit_of_interest)
    data_pre_trans <- data
    if (!is.null(data_transformation)) {
        print("Transforming data...")
        if (data_transformation == "zscore") {
            data <- scale(data)
            if (any(is.na(data))) {
                print("NA values found in data because of no variance/ all zeroes, replacing with 0...")
                data[is.na(data)] <- 0
            }
            data[data > 3] <- 3
            data[data < -3] <- -3

        } else if (data_transformation == "log10") {
            data <- log10(data + pseudocount)
        } else {
            print("Unknown data transformation, exiting function")
            return()
        }
    }

    data <- data %>% as.data.frame()
    data_metrics <- data.frame(matrix(ncol = 0, nrow = dim(data)[1]))
    data_metrics$sum_copyn <- apply(data_pre_trans, 1, \(x) sum(x))
    data_metrics$sum_uniq_f <- apply(data_pre_trans, 1, \(x) sum(x > 0))
    data_metrics$var <- apply(data, 1, var)
    # data_metrics$mean_dist <- apply(as.matrix(dist(data, method = "euclidean")), 1, mean, na.rm = TRUE)
    dists <- as.matrix(dist(data, method = "euclidean"))
    diag(dists) <- NA
    data_metrics$mean_dist <- apply(dists, 1, mean, na.rm = TRUE)
    data_metrics$`log10_#_genomes` <- data_all$`log10_#_genomes`
    # data_metrics$mOTU_ID <- data_all$mOTU_ID
    data_metrics$Family <- data_all$Family
    row_annotation <- data_metrics

    # plot data as heatmap using pheatmap package.
    print("Plotting heatmap...")
    col_annotation <- cazyAnnots %>%
        rename(cazy_family = Subfamily) %>%
        inner_join(data.frame(cazy_family = colnames(data)))
    col_annotation <- col_annotation[match(colnames(data), col_annotation$cazy_family), ]
    col_annotation <- col_annotation %>%
        as.data.frame()
    rownames(col_annotation) <- NULL
    col_annotation <- col_annotation %>%
        column_to_rownames("cazy_family")
    stopifnot(all(colnames(data) == col_annotation$cazy_family))
    col_annotation <- col_annotation %>%
        mutate(across(colnames(.), \(x) factor(x, c("Yes", "No", "NA"))))
    annot_colors_col <- map(col_annotation, \(x) return(c(Yes = 'black', No = "grey", "NA" = 'white'))) # Weirdly, does not work with colnames(col_annotation) :)
    annot_legend_helper <- map(col_annotation, \(x) return(FALSE))
    annot_colors <- annot_colors_col
    if ("Family" %in% colnames(row_annotation)) {
        row_annotation$Family <- map_chr(row_annotation$Family, \(x) str_c(str_split(x, " ")[[1]][1:2], collapse = " "))
    }
    set.seed(123123231) # for consistency of the colors in the row/col annotations
    # if (unit_of_interest == "mOTU_ID") {
    #     row_annotation <- row_annotation
    # } else {
    #     row_annotation <- NA
    # }
    # browser()
    # print(dim(data))
    # pheatmap_object <- pheatmap(
    #     data,
    #     cluster_rows = T,
    #     cluster_cols = T,
    #     show_rownames = unit_of_interest == "mOTU_ID",
    #     show_colnames = T,
    #     annotation_col = col_annotation,
    #     # annotation_row = ifelse("number_genomes" %in% colnames(data), row_annotation, NA),
    #     annotation_row = row_annotation,
    #     annotation_colors = annot_colors,
    #     # annotation_legend = F,
    #     # annotation_legend = list(DF = FALSE),
    #     fontsize = 6,
    #     breaks = seq(-3, 3, length.out = 101))
    # col_annotation <- HeatmapAnnotation(df = col_annotation, col = list(field = c("white", "black")))
    col_annotation <- HeatmapAnnotation(df = col_annotation, col = list(
        DF = c("No" = "grey", "Yes" = "black", "NA" = "white"),
        Glycogen = c("No" = "grey", "Yes" = "black", "NA" = "white"),
        PG = c("No" = "grey", "Yes" = "black", "NA" = "white"),
        GAG = c("No" = "grey", "Yes" = "black", "NA" = "white"),
        Mucin = c("No" = "grey", "Yes" = "black", "NA" = "white")
    ))
    library(circlize)
    testttttt <- map(row_annotation, \(x) {
        if (is.numeric(x)) {
            colorRamp2(c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)), c('white', 'black'))
        } else {
            NULL
        }
    })
    testttttt <- testttttt[map_lgl(testttttt, \(x) !is.null(x))]
    row_annotation <- rowAnnotation(df = row_annotation, col = testttttt)
    if (data_transformation == 'zscore') {
        color_mapping <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
    } else {
        color_mapping <- colorRamp2(c(0, max(data)), c("white", "red"))
    }


    pheatmap_object <- Heatmap(
        data,
        top_annotation = col_annotation,
        left_annotation = row_annotation,
        col = color_mapping,
        heatmap_legend_param = list(
            at = if (data_transformation == "zscore") {pretty(c(-3, 0, 3))} else {pretty(c(0, max(data)))}
        )
    )
    # color_bar = "continuous"))
    return(pheatmap_object)
}
