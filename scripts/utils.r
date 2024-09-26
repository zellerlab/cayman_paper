library(tidyverse)
library(ape)
library(RColorBrewer)
library(compositions)
library(ggrepel)
library(circlize)

cazyFamilyOrder <- c("GT", "CE", "PL", "GH", "CBM")
familyColors <- c("#47484C", #GT 
                  "#806F6F", #CE
                  "#655399", #PL
                  "#478E87", #GH
                  "#B55B9E") #CBM

mucin_pathway_colors <- c(
  "Bacteroides (7293)" = "#4a3300",
  "Barnesiella (1236)" = "#7a5a10",
  "Coprobacter (119)" = "#9e7923",
  "Parabacteroides (2953)" = "#c79f42",
  "Paraprevotella (515)" = "#f5cc6e",
  "Eisenbergiella (53)" = "#008a29",
  "Hungatella (66)" = "#2fba58",
  "Akkermansia (1313)" = "#FF61C3"
)

qq <- function(d, co, pos = 1) {
    return(d %>%
        pull({{ co }}) %>% magrittr::extract2(pos))
}

get_family_sharing_rate <- function(df, taxa, lev = "Genus"){

  jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    u = length(a) + length(b) - intersection
    return (intersection/u)
  }

  # Subset to only genera of interest.
  df <- df[df[[lev]] %in% taxa, ]

  # count ALL FAMILIES (and then later subtract)
  FamilyNs <- table(unlist(df$presentFamilies))

  specificFamilies <- list()
  for (g in unique(df[[lev]])){
    familiesInTaxon <- df %>%
      filter(.data[[lev]] == g) %>%
      pull(presentFamilies) %>%
      magrittr::extract2(1)
    tmp <- list()
    for (g2 in unique(df[[lev]])){
      familiesInTaxon2 <- df %>%
        filter(.data[[lev]] == g2) %>%
        pull(presentFamilies) %>%
        magrittr::extract2(1)
      # Calculate jaccard distance
      tmp[[length(tmp) + 1]] <- jaccard(familiesInTaxon, familiesInTaxon2)
      names(tmp)[length(tmp)] <- str_c(g, g2, sep = "__", collapse = "__")
    }
    specificFamilies[[length(specificFamilies) + 1]] <- tmp
  }
  return(specificFamilies)
}

re <- function(x) {
  # Remove PC
  #return(10^x - 1E-5)
  # Keep PC...
  return((10^x))
}

get_pw_assocs <- function(matrix1, 
                          matrix2, 
                          generaToTest = mods$genusRaw, 
                          cazyFamiliesToTest = NULL) {
  resLM <- list()
  #for (genus in colnames(pPredGenus)[!colnames(pPredGenus) %in% c("sampleID", "study_name", "non_westernized", "Unique_CAZyme_count")]) {
  for (genus in generaToTest) {
    print(genus)
    # show only genera associated with richness
    genus <- str_replace_all(genus, '`', '')
    for(cazyFamily in cazyFamiliesToTest) {
      # print
      g <- cbind(matrix1[, genus, drop = F], t(matrix2[cazyFamily, , drop = F]))
      #colnames(g)[1] <- str_replace_all(g, " ", "__")
      #print(head(g))
      print(str_c(cazyFamily, " ~ ", str_c('`', genus, "`", sep = "", collapse = "")))
      #tmp <- lm(data = g, formula = str_c(cazyFamily, " ~ ", str_c('`', genus, "`", sep = "", collapse = "")))

      tmp <- lm(data = g, formula = str_c(cazyFamily, " ~ ", genus))
      resLM[[length(resLM) + 1]] <- tmp
      names(resLM)[length(resLM)] <- str_c(genus, cazyFamily, sep = "__", collapse = "__")
    }
  }
  

  # all
  # for(cazyFamily in cazyFamiliesToTest) {
  #   g <- cbind(matrix1[, generaToTest, drop = F], t(matrix2[cazyFamily, , drop = F]))
  #     tmp <- lm(data = g, formula = str_c(cazyFamily, " ~ ", str_c(generaToTest, collapse = " + ")))
  #     resLM[[length(resLM) + 1]] <- tmp
  #     names(resLM)[length(resLM)] <- str_c("all_genera", cazyFamily, sep = "__", collapse = "__")
  # }
  
  resLM <- tibble(l = resLM,
                  namesRaw = names(resLM)) %>%
    mutate(ls = map(l, summary)) %>%
      mutate(pval = map_dbl(ls, function(s) s$coefficients[rownames(s$coefficients) != "(Intercept)", ][4])) %>%
      filter(!is.na(pval)) %>%
      mutate(beta = map_dbl(ls, function(s) s$coefficients[rownames(s$coefficients) != "(Intercept)", ][1])) %>%
      #mutate(pval = ifelse(str_detect(cazyFamily, "all_genera"), NA, pval)) %>%
      #mutate(beta = ifelse(str_detect(cazyFamily, "all_genera"), NA, beta)) %>%
    mutate(genus = str_split_fixed(namesRaw, "__", n = 2)[, 1]) %>%
    mutate(family = str_split_fixed(namesRaw, "__", n = 2)[, 2]) %>%
    arrange()
  
  return(resLM)
}


get_pw_assocs_species_level_per_genus <- function(listOfMatrices, 
                          matrix2, 
                          generaToTest = mods$genusRaw, 
                          cazyFamiliesToTest = NULL) {
  resLM <- list()
  for (genus in generaToTest) {
    # show only genera associated with richness
    genus <- str_replace_all(genus, '`', '')
    for(cazyFamily in cazyFamiliesToTest) {
      g <- cbind(listOfMatrices[[genus]], t(matrix2[cazyFamily, , drop = F]))
      print(str_c(cazyFamily, " ~ ", str_c(colnames(listOfMatrices[[genus]]), collapse = " + ", sep = " + ")))
      if (length(colnames(listOfMatrices[[genus]])) == 0) {
        next
      }
      tmp <- lm(data = g, formula = str_c(cazyFamily, " ~ ", str_c(colnames(listOfMatrices[[genus]]), collapse = " + ", sep = " + ")))
      resLM[[length(resLM) + 1]] <- tmp
      names(resLM)[length(resLM)] <- str_c(genus, cazyFamily, sep = "__", collapse = "__")
    }
  }

  
  resLM <- tibble(l = resLM,
                  namesRaw = names(resLM)) %>%
    mutate(ls = map(l, summary)) %>%
      #mutate(pval = map_dbl(ls, function(s) s$coefficients[rownames(s$coefficients) != "(Intercept)", ][4])) %>%
      mutate(pval = map_dbl(ls, function(s) min(s$coefficients[rownames(s$coefficients) != "(Intercept)", 4]))) %>%
      #filter(!is.na(pval)) %>%
      #mutate(beta = map_dbl(ls, function(s) s$coefficients[rownames(s$coefficients) != "(Intercept)", ][1])) %>%
      mutate(beta = map_dbl(ls, function(s) mean(s$coefficients[rownames(s$coefficients) != "(Intercept)", 1] > 0))) %>%
      #mutate(pval = ifelse(str_detect(cazyFamily, "all_genera"), NA, pval)) %>%
      #mutate(beta = ifelse(str_detect(cazyFamily, "all_genera"), NA, beta)) %>%
    mutate(genus = str_split_fixed(namesRaw, "__", n = 2)[, 1]) %>%
    mutate(family = str_split_fixed(namesRaw, "__", n = 2)[, 2]) %>%
    arrange()
  
  return(resLM)
}

CalcGFC <- function(x, conditionCol = 'non_westernized', levs = c("Non_Westernized", "Westernized"), probs.fc = seq(.05, .95, .05)) {
    x.pos <- x$relAb[x[[conditionCol]] == levs[1]]
    x.neg <- x$relAb[x[[conditionCol]] == levs[2]]
    q.p <- quantile(x.pos, probs = probs.fc)
    q.n <- quantile(x.neg, probs = probs.fc)
    return(sum(q.p - q.n) / length(q.p))
}

get_scatter <- function(g, f, num) {

  tmp <- test2ClrOnCountsBySampleRPKMlog10_NW %>% filter(genus == g, family == f) %>% pull(ls) %>% magrittr::extract(1)
  tmp <- tmp[[1]]
  tmp2 <- test2ClrOnCountsBySampleRPKMlog10_W %>% filter(genus == g, family == f) %>% pull(ls) %>% magrittr::extract(1)
  tmp2 <- tmp2[[1]]

  res <- get_western_scatter(
    matrix1 = pPredGenusCountsWithoutPseudoCountCLR, 
    matrix2 = nwclog10,
    genus = g,
    cazyFamily = f, 
    metaFile = meta
  )

  print(res[[2]])
  print(res[[3]])
  print(res[[4]])
  print(res[[5]])

  xRAN <- res[[2]] - res[[4]]
  yRAN <- res[[5]] - res[[3]]

  

  # cell_fun_W = colorRamp2(c(0, max(basee[cell_group == "W_Higher" & !is.na(cell_group == "W_Higher")], na.rm = TRUE)), c("white", '#3B6FB6'))
  # cell_fun_NW = colorRamp2(c(0, max(basee[cell_group == "NW_Higher" & !is.na(cell_group == "NW_Higher")], na.rm = TRUE)), c("white", '#D41645'))

  # print(length(res))
  return(res[[1]] +
    #geom_abline(intercept = tmp$coefficients[1, 1], slope = tmp$coefficients[2, 1], color = '#D41645') +
    geom_abline(intercept = tmp$coefficients[1, 1], slope = tmp$coefficients[2, 1], color = '#D41645') +
    annotate(geom = 'text', x = res[[4]] + 0.5 * xRAN + -0.03 * xRAN, y = res[[5]] - 0.99 * yRAN, label = bquote("R"^"2" * " (NW): " ~ .(round(tmp$r.squared, 2))), color = '#D41645', hjust = 0, size = 3.25) +
    #geom_abline(intercept = tmp2$coefficients[1, 1], slope = tmp2$coefficients[2, 1], color = '#3B6FB6') +
    geom_abline(intercept = tmp2$coefficients[1, 1], slope = tmp2$coefficients[2, 1], color = '#3B6FB6') +
    annotate(geom = 'text', x = res[[4]] + 0.5 * xRAN + -0.03 * xRAN, y = res[[5]] - 0.9 * yRAN, label = bquote("R"^"2" * " (W): " ~ .(round(tmp2$r.squared, 2))), color = '#3B6FB6', hjust = 0, size = 3.25) +
    annotate(geom = 'text', x = res[[4]] + -0.01 * xRAN, y = res[[5]] - 0.025 * yRAN, label = num, color = 'black', hjust = 0, size = 5.5)
  )
}

get_western_scatter <- function(matrix1,
                                matrix2,
                                genus = mods$genusRaw,
                                cazyFamily = NULL,
                                metaFile = NULL) {
  resLM <- list()
  genus <- str_replace_all(genus, '`', '')
  g <- cbind(matrix1[, c('sampleID', genus), drop = F], t(matrix2[cazyFamily, , drop = F]))
  g <- inner_join(g, metaFile %>% select(sampleID, non_westernized), by = 'sampleID')
    if (cazyFamily == "GH95") {
    # Family has 2 outliers with log10 < 1, which make te Figure less informative.
    # I'm removing them here
    g <- g %>%
        filter(GH95 > 1) 
  }
  colnames(g) <- c("sampleID", "genus", "cazyFamily", "westernization")
  g <- g %>%
    mutate(al = ifelse(westernization == "Non_Westernized", 0.05, 0.01))
  print(g %>% head())
  print(min(g[['cazyFamily']]))
  g <- g %>%
    rename(westernized = westernization) %>%
    mutate(westernized = ifelse(westernized == "Non_Westernized", "No", "Yes")) %>%
    mutate(westernized = factor(westernized, levels = c("No", "Yes")))
#cell_fun_W = colorRamp2(c(0, max(basee[cell_group == "W_Higher" & !is.na(cell_group == "W_Higher")], na.rm = TRUE)), c("white", '#3B6FB6'))
#cell_fun_NW = colorRamp2(c(0, max(basee[cell_group == "NW_Higher" & !is.na(cell_group == "NW_Higher")], na.rm = TRUE)), c("white", '#D41645'))    
  p <- ggplot(g, aes(x = genus, y = cazyFamily, color = westernized, alpha = al)) +
    # p <- ggplot(g, aes(x = genus, y = cazyFamily, alpha = al)) +
    geom_point() +
    theme_classic() +
    #ggtitle(str_c(genus, '\n', cazyFamily, sep = "", collapse = "")) +
    scale_alpha_continuous(range = c(0.05, 0.1)) +
    scale_color_manual(values = c("Yes" = "#3B6FB6", "No" = "#D41645")) +
    xlab(str_c(genus, " (CLR)")) +
    ylab(str_c(cazyFamily, " (log10)")) +
    guides(alpha = "none", color = "none")
  
  return(list(p, max(g[['genus']]), min(g[['cazyFamily']]), min(g[['genus']]), max(g[['cazyFamily']])))
}


overlapCounterBothCardinalities <- function(queryFamily, ref = NULL) {
  # totalNumberOfSequencesThatContainsFamily queryFamily
  #n <- totalCAZymeCounts %>% filter(cazy_family == queryFamily) %>% pull(n)
  # Sort co-occuring families
  if (!queryFamily %in% colnames(ref)) {
    return(NA)
  }
  tmp <- ref
  # tmp <- tmp %>%
  #   column_to_rownames('genome') %>%
  #   as.data.frame()
  co <- tmp[[queryFamily]] == 1
  r <- list()
  for (family in colnames(tmp)) {
    coDown <- tmp[[family]] == 1
    # Normalize by the cardinality of the 'queryFamily'
    r[[length(r) + 1]] <- sum(co & coDown) / sum(co)
  }
  names(r) <- colnames(tmp)
  return(unlist(r))
}

get_genera_assocs_with_crc_cbms <- function(baseData = NULL, taxonName = NULL) {
  genera_assocs_with_crc_cbms <- list()
  taxa <- list()
  for (i in 1:length(baseData[[taxonName]])) {  
    genus <- baseData[[taxonName]][i]
    data <- baseData$oFin[[i]]
    genera_assocs_with_crc_cbms[[length(genera_assocs_with_crc_cbms) + 1]] <- rbind(
      inner_join(data, CRC_assoc_cazymes, by = c('inner' = 'cazyme_family')),
      inner_join(data, CRC_assoc_cazymes, by = c('outer' = 'cazyme_family'))
    ) %>%
      # This is not a good idea anymore...
      #filter(jaccard > 0.2) %>%
        filter(totalN.inner > 5 & totalN.outer > 5)
    genera_assocs_with_crc_cbms[[length(genera_assocs_with_crc_cbms)]][[taxonName]] <- genus
  }

  plotData <- do.call('rbind', genera_assocs_with_crc_cbms)
  #plotData[[taxonName]] <- taxa
  plotData <- plotData %>%
    distinct(inner, outer, jaccard, Genus, .keep_all = TRUE) %>%
    mutate(pairs = map2_chr(inner, outer, function(i, o) {
      return(str_c(sort(c(i, o)), sep = "__", collapse = "__"))
    })) %>%
    group_by(Genus, pairs) %>%
      nest() %>%
      mutate(data = map(data, function(x) {
        if (any(x$jaccard > 0.2)) {
          return(x %>% arrange(desc(jaccard)) %>% head(1))
        } else {
          return(NULL)
        }
      })) %>%
      filter(!is.null(data)) %>%
      unnest() %>%
    ## This used to work (when we only normalized wrt to one cardinality)
    ## THIS DISTINCT CALL INTRODUCED A MAJOR BUG!
    #distinct(pairs, .keep_all = TRUE) %>%
    # There is one cbm-cbm pair that messes the plot up.
    # I think we'll ultimately not show it but for now:
    #mutate(pairs = ifelse(pairs == "CBM32__CBM51", "CBM51__CBM32", pairs)) %>%
    # Fuck it, those are not interesting
    filter(! (str_detect(inner, "CBM") & str_detect(outer, "CBM"))) %>%
    mutate(innerPlot = str_split_fixed(pairs, "__", n = 2)[, 1]) %>%
    mutate(outerPlot = str_split_fixed(pairs, "__", n = 2)[, 2]) %>%
    mutate(N = (totalN.inner + totalN.outer) / 2) %>%
    filter(!str_detect(.data[[taxonName]], "^NA"))

  plotData <- plotData %>%
    left_join(cazyAnnots, by = c('outer' = "Subfamily")) %>%
    mutate(substrateClass =
      case_when(
        (str_detect(FUNCTION_AT_DESTINATION_1, "Mucin") | str_detect(FUNCTION_AT_DESTINATION_1, "mucin")) & (str_detect(FUNCTION_AT_DESTINATION_1, "GAG") | str_detect(FUNCTION_AT_DESTINATION_1, "gag")) ~ "both",
        str_detect(FUNCTION_AT_DESTINATION_1, "Mucin") | str_detect(FUNCTION_AT_DESTINATION_1, "mucin") ~ "Mucin",
        str_detect(FUNCTION_AT_DESTINATION_1, "GAG") | str_detect(FUNCTION_AT_DESTINATION_1, "gag") ~ "GAG",
        .default = NA))

  library(ggunileg)

  p <- ggplot() +
    geom_point(data = plotData, aes_string(x = "innerPlot", y = "jaccard", size = "N", color = taxonName, shape = taxonName, group = 'pairs'), position = 'dodge') +
    geom_point(data = plotData %>% filter(!is.na(substrateClass)), aes_string(x = "innerPlot", y = "jaccard"), size = 5, shape = 1, color = 'red') +
    geom_text_repel(data = plotData, aes(x = innerPlot, y = jaccard, label = outerPlot), max.overlaps = Inf) +
    scale_color_highres(num_shape_levels = 5, name = taxonName) +
    guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) + 
    theme_classic() +
    xlab("CBMs") +
    ylab("Degree of co-occurence\n[within ORF]")
  #ggsave(plot = p, filename = outPath, width = 12, height = 5)
  return(plotData)

}
