library(here)
source(here('scripts', 'utils.r'))
library(stringr)
library(igraph)
library(readxl)


#crc_assoc_motus <- readRDS(here('data', 'Nic_taxonomy_CRC.RDS')) %>%
crc_assoc_motus <- read_csv(here('data', 'Intermediate_Files', 'taxonomy_siamcat_CRC.csv')) %>%
as_tibble() %>%
filter(fc > 0) %>%
filter(p.adj < 0.01) %>% arrange(desc(fc)) %>%
#print(n=25) %>%
head(25) %>%
select(mOTU_nr) %>%
rename(mOTU_ID = mOTU_nr)

# Load cazy annotations
#cazyAnnots <- read_tsv(here('data', "20230609_glycan_annotations_cleaned_Nic.tsv"))
completed_substrate_annotations <- read_xlsx(here("data", "Glycan_Annotations", "20230607_glycan_annotations_cleaned.xlsx"))
glycan_annotations_final_cleaned <- completed_substrate_annotations %>% select(c(Family:Subfamily,ORIGIN:FUNCTION_AT_DESTINATION_3, Glycan_annotation))
cazyAnnots <- glycan_annotations_final_cleaned

#cazy_richness <- read_tsv(here('data', "20230929_CAZyme_count_df_Nic.tsv")) %>%
cazy_richness <- read_tsv(here('data', "Intermediate_Files", "20240809_CAZyme_count_df.tsv")) %>%  
  rename(sampleID = Sample) %>%
  select(all_of(c("sampleID", "Unique_CAZyme_count")))

motus3.0_taxonomy <- read_tsv(here('data', "motus3.0_taxonomy.tsv"))

almeidaCAZy <- read_tsv(here('data', 'almeida_cazy_annotations.tsv'))

tree.filtered <- read.tree(here('data', "tree.genus.ncbi.filtered.nwk"))

CRC_assoc_cazymes <- data.frame(cazyme_family = c("CBM66","CBM51", "CBM62", "CBM16", "CBM70", "CBM40"))

onlySeqsWithMultipleDifferentFamilies <- almeidaCAZy %>%
  group_by(sequenceID, genome) %>%
  filter(n() > 1)

onlySeqsWithMultipleDifferentFamilies <- onlySeqsWithMultipleDifferentFamilies %>%
  filter(length(unique(cazy_family)) > 1)


d <- function(x) {
  return(ifelse(sum(x) >= 1, 1 , 0))
}

onlySeqsWithMultipleDifferentFamiliesFinal <- onlySeqsWithMultipleDifferentFamilies %>%
  # unnest() %>%
  mutate(bla = 1) %>%
  pivot_wider(id_cols = c(sequenceID, Lineage, Kingdom, Phylum, Class, Order, Family, Genus, Species, mOTU_ID), names_from = cazy_family, values_from = bla, values_fill = 0, values_fn = d)

onlySeqsWithMultipleDifferentFamiliesFinalByGenus <- onlySeqsWithMultipleDifferentFamiliesFinal %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>%
  nest() %>%
  # filter(
  #   str_detect(Species, "Fuso") |
  #   str_detect(Species, "nucleatum") |
  #     str_detect(Species, "Parvimonas") |
  #     str_detect(Species, "Peptostreptococcus") |
  #     str_detect(Species, "Porphyromonas")) %>%
  # filter(!str_detect(Species, "sedis")) %>%
  left_join(almeidaCAZy %>% select(genome, Genus) %>% distinct() %>% group_by(Genus) %>% tally() %>% rename(numberMAGs = n), by = 'Genus') %>%
  filter(numberMAGs > 10)

ofInterestGenusBothCardinalities <- onlySeqsWithMultipleDifferentFamiliesFinalByGenus %>%
#ofInterestBySpecies <- onlySeqsWithMultipleDifferentFamiliesFinalBySpecies %>%
#ofInterestBySpeciesCBMsAsReferenceSet <- onlySeqsWithMultipleDifferentFamiliesFinalBySpecies %>%  
  #filter(str_detect(Genus, "Fusobacterium")) %>%
  #head(1) %>%
  # filter(str_detect(Genus, "Fusobacterium") | 
  # str_detect(Genus, "Porphyromonas") |
  # str_detect(Genus, "Peptostreptococcus") |
  # str_detect(Genus, "Parvimonas")) %>% 
  mutate(o = map(data, function(x) {
    print(dim(x))
    r <- list()
    u <- colnames(x)[colnames(x) != 'sequenceID' & colnames(x) != 'Lineage']
    #u <- colnames(x)[!colnames(x) %in% c("sequenceID", "Lineage", "Species", "mOTU_ID")]
    for (family in u) {
      r[[length(r) + 1]] <- overlapCounterBothCardinalities(family, x)
    }
    names(r) <- u
    rD <- map(r, as.data.frame)
    rD <- map2(rD, names(r), function(x, y) {
      # outer should be here the family with the reference cardinality
      return(x %>% mutate(outer = y))
    })
    finalCoOccurancesWithinORFs <- do.call('rbind', rD) %>%
      rownames_to_column('inner') %>%
      as_tibble() %>%
      rename(jaccard = `.x[[i]]`) %>%
      mutate(inner = stringr::str_split_fixed(inner, "[.]", n = 2)[, 2]) %>%
      filter(outer != inner) %>%
      #inner_join(to, by = c('inner' = "family")) %>%
      #inner_join(to, by = c('outer' = "family"), suffix = c(".inner", ".outer"))
      arrange(desc(jaccard))
    return(finalCoOccurancesWithinORFs)
  }))

ofInterestGenusBothCardinalities <- ofInterestGenusBothCardinalities %>%
  mutate(oFin = map2(o, data, \(x, d) {
    to <- d %>%
      select(-any_of(c('sequenceID', 'Lineage', 'Species', 'mOTU_ID'))) %>%
      apply(., 2, sum) %>%
      as.data.frame() %>%
      rownames_to_column('family') %>%
      rename(totalN = '.')
    return(x %>%
      arrange(desc(jaccard)) %>%
      #filter(jaccard > 0.05) %>%
      inner_join(to, by = c('inner' = "family")) %>%
      inner_join(to, by = c('outer' = "family"), suffix = c(".inner", ".outer"))
      ) 
  }))

tmp2 <- get_genera_assocs_with_crc_cbms(ofInterestGenusBothCardinalities, "Genus")
tmp2

library(network)
library(sna)
library(ggplot2)
library(RColorBrewer)

edgeWeightFormula <- function(bla) {
  o <- sqrt((((bla))))
  return(o)
}

set.seed(1233)
vertexColors <- c(sample(brewer.pal(4, "Set1"), 4), "grey")
vertexLabels <- c("Mucin,GAG", "Mucin", "GAG", "Fibre", "Other")

substrateColors <- c(sample(brewer.pal(4, "Set1"), 4), "grey")
names(substrateColors) <- c("Mucin,GAG", "Mucin", "GAG", "Fibre", "Other")

onlyCRCGenera <- FALSE
hmBool <- FALSE

networks <- tmp2 %>%
# TODO: Filter by CRC Genera
mutate(Genus = str_split_fixed(Genus, " ", n = 2)[, 2]) %>%
{if (onlyCRCGenera) {
(.) %>%
inner_join(
    #read_csv('/g/scb2/zeller/pekel/CRC_meta_idtaxa/Results/Association_table_ncbi_16SLOCRC.csv') %>% 
    #read_csv('/g/scb2/zeller/karcher/CAZY_project_v2/data/Association_table_ncbi_all_data.csv') %>%
    read_csv(here('data', 'crc_enriched_genera.tsv')) %>%
    rename(genus = `...1`) %>% 
    filter(fc > 0) %>% 
    filter(p.adj < 0.1) %>% 
    select(genus) %>% 
rename(Genus = genus))
} else {
(.)
}} %>%
# innerPlot corresponds to the CBM, i.e. central node
group_by(innerPlot) %>%
nest() %>%
#filter(innerPlot == "CBM66") %>%
mutate(dataNetworkRaw = map2(data, innerPlot, \(x, i) {
x <- x %>%
    ungroup() %>%
    #mutate(Genus = str_split_fixed(Genus, " ", n = 2)[, 2]) %>%
    mutate(genusOuterPlot = str_c(Genus, "___", outerPlot)) %>%
    select(Genus, genusOuterPlot, jaccard)
levs <- unique(c(i, x$genusOuterPlot))
levs2 <- map_chr(levs[2:length(levs)], \(x) str_split(x, "___")[[1]][2])
levs <- c(levs[str_detect(levs, "CBM")], levs[2:length(levs)][order(levs2)])
    x <- x %>%
    mutate(i = factor(i, levels = levs)) %>%
    mutate(genusOuterPlot = factor(genusOuterPlot, levels = levels(i))) %>%
    pivot_wider(id_cols = genusOuterPlot, names_from = i, values_from = jaccard, names_expand = TRUE) %>%
    column_to_rownames('genusOuterPlot')
stopifnot(str_detect(colnames(x)[1], "CBM"))
x <- rbind(rep(NA, (dim(.)[2])), x)
rownames(x)[1] <- i
x <- x[, levs]
x <- x[levs, ]
stopifnot(all(rownames(x) == colnames(x)))
return(x)
}))  %>%
mutate(adjacencyGraphQuant = map(dataNetworkRaw, \(x) {
x <- as.matrix(x)
x[is.na(x)] <- 0
return(x)
})) %>%
mutate(adjacencyGraphQuantOrigData = map(dataNetworkRaw, \(x) {
x <- as.matrix(x)
x[is.na(x)] <- 0
return(x)
})) %>%    
# For bipartite
mutate(adjacencyGraphQuant = map(dataNetworkRaw, \(x) {
x <- as.matrix(x)
x[is.na(x)] <- 0
rownames(x) <- c(rownames(x)[1], str_split_fixed(rownames(x)[2:length(rownames(x))], "___", n = 2)[, 2])
colnames(x) <- c(colnames(x)[1], str_split_fixed(colnames(x)[2:length(colnames(x))], "___", n = 2)[, 2])
print(x)
x <- x[!duplicated(rownames(x)), ]
x <- x[, !duplicated(colnames(x)) ]
stopifnot(all(rownames(x) == colnames(x)))
return(x)
})) %>%
mutate(finalNetworks = map(adjacencyGraphQuant, \(a) {
bla <- matrix(0, nrow = dim(a)[1], ncol = dim(a)[2])
bla[, 1] <- 1
g <- graph_from_adjacency_matrix(
    bla,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
)
# V(g)$cazyName <- map_chr(rownames(a), \(x) str_split(x, "___")[[1]][2])
V(g)$cazyName <- rownames(a)
V(g)$cazyName[1] <- rownames(a)[1]
#E(g)$taxon <- c(map_chr(rownames(a)[2:length(rownames(a))], \(x) str_split(x, "___")[[1]][1]))
#stopifnot(length(a[, 1]) - 1 == length(edge_attr(g)$weight))
#E(g)$weight <- a[2:(dim(a)[1]), 1]
V(g)$type <- ifelse(str_detect(V(g)$cazyName, "CBM"), TRUE, FALSE)
return(g)
})) %>%
mutate(layouts = map2(finalNetworks, adjacencyGraphQuant, \(g, a) {
# LO = layout_as_star(g)
LO = layout_as_bipartite(g)
#LO = (LO * ((1 - edgeWeightFormula(a[, 1])) + 1))
# LO = LO * sample(1:5, size = length(a[, 1]), replace = TRUE)
return(LO)
# })) %>%
#   mutate(edgeColorMap = map())
})) %>%
mutate(nodeColorMap = map(adjacencyGraphQuant, \(a) {
data <- a[, 1]
# data <- map_chr(names(data)[2:length(names(data))], \(x) str_split(x, "___")[[1]][2])
data <- names(data)
data <- unique(data)
data <- data.frame(cazy_family = data) %>%
    left_join(cazyAnnots, by = c('cazy_family' = "Subfamily")) %>%
    mutate(substrateClass =
    case_when(
        (str_detect(FUNCTION_AT_DESTINATION_1, "Mucin") | str_detect(FUNCTION_AT_DESTINATION_1, "mucin")) & (str_detect(FUNCTION_AT_DESTINATION_1, "GAG") | str_detect(FUNCTION_AT_DESTINATION_1, "gag")) ~ "Mucin,GAG",
        str_detect(FUNCTION_AT_DESTINATION_1, "Mucin") | str_detect(FUNCTION_AT_DESTINATION_1, "mucin") ~ "Mucin",
        str_detect(FUNCTION_AT_DESTINATION_1, "GAG") | str_detect(FUNCTION_AT_DESTINATION_1, "gag") ~ "GAG",
        .default = "Other")) %>%
    select(cazy_family, substrateClass) %>%
    mutate(substrateClass = factor(substrateClass, levels = vertexLabels)) %>%
    mutate(cols = map_chr(as.numeric(substrateClass), \(x) vertexColors[x])) %>%
    select(cazy_family, cols)
return(data)
})) %>%
mutate(nodeColors = map2(finalNetworks, nodeColorMap, function(g, ncm) {
return(data.frame(cazy_family = vertex_attr(g)$cazyName) %>%
    left_join(ncm, by = 'cazy_family') %>%
    pull(cols))
})) %>%
mutate(heatmapBase = map(adjacencyGraphQuantOrigData, \(x) {
#print(x)
x <- x[-1, , drop = F]
#print(x)
x <- x[, 1, drop = F]
#print(x)
genus <- str_split_fixed(rownames(x), "___", n = 2)[, 1]
family <- str_split_fixed(rownames(x), "___", n = 2)[, 2]
return(data.frame(genus = genus, family = family, jaccard = x[, 1]))
}))

g <- networks %>% qq(finalNetworks)
d <- networks %>% qq(adjacencyGraphQuant)


# pts.circle[, 1] <- map2_dbl(pts.circle[, 1], pts.circle[, 2], \(x, y) x * cos(te) - y * sin(te))
# pts.circle[, 2] <- map2_dbl(pts.circle[, 1], pts.circle[, 2], \(x, y) x * sin(te) + y * cos(te))

pmap(list(networks$finalNetworks, networks$innerPlot, networks$layouts, networks$nodeColors), \(g, a, l, nc) {
width <- height <- 7
# nn <- length(V(g)$cazyName)
# n <- (nn*100) # number of points you want on the unit circle
# pts.circle <- t(sapply(1:n, function(r) c(cos(2 * r * (pi) / n), sin(2 * r * (pi) / n))))
# pts.circle <- pts.circle[seq(1, nn*100, 100) + 1, ]
if (!onlyCRCGenera && a == "CBM51") {
height <- 12
width <- 12
}
# pdf(file = str_c("/g/scb2/zeller/karcher/CAZY_project_v2/plots/", a, "_co_occurence_networks_only_CRC_genera_", onlyCRCGenera, ".pdf"), width = width, height = height)
# bla <- plot(
#   g,
#   layout = l,
#   # vertex.label.dist = 2.5,
#   vertex.label.color = "black",
#   # edge.color = adjustcolor('darkgrey', alpha.f = (E(g)$weight)),
#   # edge.color = map_chr(E(g)$weight, \(x) adjustcolor('darkgrey', x)),
#   vertex.color = nc,
#   # edge.width = (E(g)$weight) * 3,
#   edge.width = 5,
#   vertex.shape = c('square', rep("circle", length(V(g)$cazyName) - 1)),
#   vertex.label = V(g)$cazyName,
#   # edge.label = E(g)$taxon
#   # edge.label.y = pts.circle[, 2] * .5,
#   # edge.label.x = pts.circle[, 1]*.5
# )
# legend("bottomright", legend = vertexLabels, pch = 21,
#   col = vertexColors, pt.bg = vertexColors, pt.cex = 1, cex = 0.8, bty = "o", ncol = 1)
# dev.off()
# return(bla)
})


data <- networks %>%
select(innerPlot, heatmapBase) %>%
unnest() %>%
pivot_wider(id_cols = c(innerPlot, family), names_from = genus, values_from = jaccard, values_fill = NA) %>%
pivot_longer(-c(innerPlot, family)) %>%
rename(genus = name, jaccard = value) %>%
left_join(motus3.0_taxonomy %>%
select(Phylum, Genus) %>%
mutate(across(colnames(.), \(x) str_split_fixed(x, " ", n = 2)[, 2])) %>% 
distinct(), by = c('genus' = "Genus")) %>%
left_join(cazyAnnots, by = c('family' = 'Subfamily')) %>%
#select(innerPlot, family, FUNCTION_AT_DESTINATION_1) %>%
mutate(substrateClass =
case_when(
(str_detect(FUNCTION_AT_DESTINATION_1, "Mucin") | str_detect(FUNCTION_AT_DESTINATION_1, "mucin")) & (str_detect(FUNCTION_AT_DESTINATION_1, "GAG") | str_detect(FUNCTION_AT_DESTINATION_1, "gag")) ~ "Mucin,GAG",
str_detect(FUNCTION_AT_DESTINATION_1, "Mucin") | str_detect(FUNCTION_AT_DESTINATION_1, "mucin") ~ "Mucin",
str_detect(FUNCTION_AT_DESTINATION_1, "GAG") | str_detect(FUNCTION_AT_DESTINATION_1, "gag") ~ "GAG",
str_detect(FUNCTION_AT_DESTINATION_1, "DF") ~ "Fibre",
.default = "Other")) 

set.seed(1233)
#substrateColors <- c(sample(brewer.pal(4, "Set1"), 3), "grey")
#names(substrateColors) <- c("Mucin,GAG", "Mucin", "GAG", "Other")

data$Phylum <- as.factor(data$Phylum)
#data$genus <- as.factor(data$genus)
data$genus <- factor(as.character(data$genus), levels = c("Bifidobacterium","Varibaculum","Alistipes","Bacteroides","Paraprevotella","Prevotella","Anaeromassilibacillus","Anaerostipes","Bacillus","Blautia","Clostridium","Coprobacillus","Dorea","Eisenbergiella","Enterococcus","Erysipelatoclostridium","Eubacterium","Faecalitalea","Flavonifractor","Hungatella","Intestinimonas","Lactobacillus","Paenibacillus","Roseburia","Ruminococcus","Ruthenibacterium","Streptococcus","Turicibacter","Tyzzerella","Fusobacterium","Brevundimonas","Stenotrophomonas","Akkermansia"))
lP <- levels(data$Phylum)
lG <- levels(data$genus)
p <- ggplot(data = data, aes(x = genus, y = family)) +
{if(!hmBool) {
list(geom_point(aes(size = jaccard, color = substrateClass)),
    scale_color_manual(values = substrateColors),
    scale_size(c(0.1, 0.5)),
    scale_size_continuous(range = c(0.1, 3)))
} else {
list(geom_tile(aes(fill = jaccard)), scale_fill_continuous(low = "thistle2", high = "darkred", na.value = "white"))
}} +
facet_grid(innerPlot ~ ., scales = 'free', space = "free", switch = "y") +
    theme_classic() +
    theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank())
    #plot.margin = unit(c(0, 0, 0, 0), "cm"),
    
    #panel.spacing.x = unit(-0.3, "lines")
#scale_x_discrete(expand = c(0, 0)) +
    NULL

# p.right <- data %>%
#   select(innerPlot, family) %>%
#   distinct() %>%
#   left_join(cazyAnnots, by = c('family' = 'Subfamily')) %>%
#   select(innerPlot, family, FUNCTION_AT_DESTINATION_1) %>%
#   mutate(substrateClass =
#   case_when(
#     (str_detect(FUNCTION_AT_DESTINATION_1, "Mucin") | str_detect(FUNCTION_AT_DESTINATION_1, "mucin")) & (str_detect(FUNCTION_AT_DESTINATION_1, "GAG") | str_detect(FUNCTION_AT_DESTINATION_1, "gag")) ~ "Mucin,GAG",
#     str_detect(FUNCTION_AT_DESTINATION_1, "Mucin") | str_detect(FUNCTION_AT_DESTINATION_1, "mucin") ~ "Mucin",
#     str_detect(FUNCTION_AT_DESTINATION_1, "GAG") | str_detect(FUNCTION_AT_DESTINATION_1, "gag") ~ "GAG",
#     .default = "Other"))

# p.right <- ggplot() +
#   geom_bar(data = p.right, aes(x = 1, y = family, fill = substrateClass), stat = 'identity') +
#   theme_classic() +
#   facet_grid(innerPlot ~ ., scales = 'free', space = "free", switch = 'y') +
#   theme(
#     axis.text.x = element_blank(),
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     plot.margin = unit(c(0,0,0,0), "cm")) +
#   scale_fill_manual(values = substrateColors)

p.bottom <- data %>%
    ungroup() %>%
    select(genus) %>%
    distinct() %>%
    left_join(
    #read_csv('/g/scb2/zeller/pekel/CRC_meta_idtaxa/Results/Association_table_ncbi_16SLOCRC.csv') %>%
    read_csv('/g/scb2/zeller/karcher/CAZY_project_v2/data/Association_table_ncbi_all_data.csv') %>%
        rename(genus = `...1`) %>%
        #mutate(crcAssociated = fc > 0 & p.adj < 0.1) %>%
        # filter(fc > 0) %>%
        # filter(p.adj < 0.1) %>%
        select(genus, fc, p.adj)) %>%
    # mutate(crcAssociated = ifelse(is.na(crcAssociated), FALSE, crcAssociated)) %>%
    mutate(crcAssociated = factor(case_when(
    fc > 0 & p.adj < 0.1 ~ "CRC-enriched",
    fc < 0 & p.adj < 0.1 ~ "CRC-depleted",
    p.adj >= 0.1 ~ "Neither"), levels = c("CRC-enriched", "CRC-depleted", "Neither"))) %>%
left_join(motus3.0_taxonomy %>%
select(Phylum, Genus) %>%
mutate(across(colnames(.), \(x) str_split_fixed(x, " ", n = 2)[, 2])) %>% 
distinct(), by = c('genus' = "Genus"))

p.bottom$genus <- factor(as.character(p.bottom$genus), levels = lG)

p.bottom <- ggplot() +
geom_bar(data = p.bottom, aes(x = genus, y = 1, fill = crcAssociated), stat = 'identity') +
#geom_point(data = p.bottom %>% mutate(bla = 1), aes(x = genus, y = 1, fill = crcAssociated)) +
#facet_grid( ~ Phylum, scales = 'free', space = "free") +
theme_classic() +
theme(
axis.text.y = element_blank(),
axis.title.y = element_blank(),
axis.ticks.y = element_blank(),
axis.text.x = element_text(angle = 45, hjust = 1),
#axis.ticks.x = element_blank(),
#axis.title.x = element_blank(),
strip.background = element_blank(),
strip.text = element_blank()) +
#plot.margin = unit(c(0,0,0,0), "cm"),
#panel.spacing.x = unit(0.1, "lines")) +
# panel.spacing = unit(0.25, "lines")) +
scale_fill_manual(values = c("CRC-enriched" = "#DE3163", "CRC-depleted" = "#4CBB17", "Neither" = 'lightgrey')) +
#scale_x_discrete(expand = c(0, 1)) +
    NULL


library(patchwork)
ggsave(plot = p + p.bottom + plot_layout(
    heights = c(30, 1),
    guides = 'collect'),
    #filename = str_c("/g/scb2/zeller/karcher/CAZY_project_v2/plots/crcheatmap_onlyCRCGenera", onlyCRCGenera, "__heatmap", hmBool, ".pdf"), 
    filename = here('figures', "Fig5_C.pdf"),
width = 7.5, height = 5.75)
