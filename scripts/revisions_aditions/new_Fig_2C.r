library(here)
library(ggrepel)
source(here('scripts', 'utils.r'))

# motus3.0_taxonomy <- read_tsv(here('data', "motus3.0_taxonomy.tsv"))
almeidaCAZy <- read_tsv(here('data', 'almeida_cazy_annotations.tsv'))
motus3.0_taxonomy <- read_tsv(here('data', "motus3.0_taxonomy.tsv"))
tree.filtered <- read.tree(here('data', "tree.genus.ncbi.filtered.nwk"))

species_of_interest <- almeidaCAZy %>%
  group_by(Species) %>%
  summarize(num_genomes = length(unique(genome))) %>%
  filter(num_genomes >= 10) %>%
  ungroup() %>%
  select(Species)

uniqueFamiliesRaw <- almeidaCAZy %>%
  # We're only interested in Species that we can name...
  filter(!str_detect(Species, "gen[.]")) %>%
  filter(!str_detect(Species, "incertae")) %>%
  filter(!str_detect(Species, "^NA")) %>%
  inner_join(species_of_interest) %>%
  group_by(mOTU_ID) %>%
  nest() 
uniqueFamiliesRaw <- uniqueFamiliesRaw %>%
  #mutate(Species = map_chr(Species, \(x) str_c(str_split(x, " ")[[1]][1:3], collapse = " ")))
  identity()
# 
uniqueFamiliesRaw <- uniqueFamiliesRaw %>% mutate(presentFamilies = map(data, function(x) {
  numGenomes <- length(unique(x$genome))
  return(x %>%
    select(genome, cazy_family) %>%
    distinct() %>%
    group_by(cazy_family) %>%
    tally() %>%
    mutate(n_N = n / numGenomes) %>%
    # Call a family present if it's found in at least 20% of the genomes.
    filter(n_N > 0.20) %>%
    pull(cazy_family))
}))

#cazyFamilyPairwiseJaccards <- get_family_sharing_rate(uniqueFamiliesRaw, tree.filtered$tip.label)
# This runs for around 20 minutes - sorry for bad code :heart:.
#cazyFamilyPairwiseJaccards <- get_family_sharing_rate(uniqueFamiliesRaw, uniqueFamiliesRaw$Species, lev = "Species")
cazyFamilyPairwiseJaccards <- get_family_sharing_rate(uniqueFamiliesRaw, uniqueFamiliesRaw$mOTU_ID, lev = "mOTU_ID")
#o <- cazyFamilyPairwiseJaccards
cazyFamilyPairwiseJaccards <- cazyFamilyPairwiseJaccards %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(jaccard_similarity = V1) %>%
  rownames_to_column('group') %>%
  mutate(taxon_1 = str_split_fixed(group, "__", n = 2)[, 1]) %>%
  mutate(taxon_2 = str_split_fixed(group, "__", n = 2)[, 2]) %>%
  pivot_wider(id_cols = taxon_1, names_from = taxon_2, values_from = jaccard_similarity) %>%
  as.data.frame()

pcoa <- cmdscale((as.dist(1-(cazyFamilyPairwiseJaccards %>%
                            #mutate(taxon_1 = str_replace_all(taxon_1, "[.]", " ")) %>%
                            column_to_rownames('taxon_1') %>%
                            as.matrix()))), k = 2)
colnames(pcoa) <- c("PCo 1", "PCo 2")
pcoa <- pcoa %>%
                                     as.data.frame() %>%
                                     mutate(mOTU_ID = rownames(.),
                                            #yOffset = `PCo 2` + 0.015) %>%
                                            yOffset = `PCo 2`) %>%           
                                     inner_join(motus3.0_taxonomy %>%
                                                 #mutate(Genus = str_split_fixed(Genus, " ", n = 2)[, 2]) %>%
                                                 mutate(Species = map_chr(Species, \(x) str_c(str_split(x, " ")[[1]][2:3], collapse = " "))) %>%
                                                 mutate(Species = str_replace(Species, "sp[.]", "sp ")) %>%
                                                 mutate(Genus = str_split_fixed(Genus, " ", n = 2)[, 2]) %>%
                                                 mutate(Phylum = str_split_fixed(Phylum, " ", n = 2)[, 2]) %>%
                                                 select(mOTUs_ID, Species, Genus, Phylum) %>% distinct() %>%
                                                 rename(mOTU_ID = mOTUs_ID,
                                                        genus = Genus,
                                                        species = Species,
                                                        phylum = Phylum), by = 'mOTU_ID') %>%
                                    group_by(phylum) %>%
                                    filter(n() > 2) %>%
                                     # Collinsella has an ill-defined phylum. Fix.
                                     #mutate(phylum = ifelse(genus == "Collinsella", "Actinobacteria", phylum)) %>%
                                     mutate(`PCo 2` = -1 * `PCo 2`) %>%
                                     mutate(`PCo 1` = 1 * `PCo 1`)
ggsave(plot = ggplot() + geom_point(data = pcoa, aes(x = `PCo 1`, y = `PCo 2`, color = phylum), alpha = 0.5) +
  #geom_text(aes(label = genus, y = yOffset)) +
  geom_text_repel(data = pcoa %>%
                    mutate(kickBool = pmap_lgl(list(`PCo 1`, `PCo 2`, species), function(po1, po2, g) {
                      # These are equivalent to the top 8 mucin-targetting taxa from 2B
                      if (g %in% c("Akkermansia muciniphila",
                                   #"Bacteroides thetaiotaomicron",
                                   #"Bacteroides vulgatus",
                                   "Bacteroides uniformis",
                                   #"Bacteroides fragilis",
                                   "Eisenbergiella tayi",
                                   "Hungatella hathewayi"
                                   #"Dialister succinatiphilus",
                                   #"Veillonella atypica",
                                   #"Eggerthella lenta",
                                   #"Bifidobacterium longum",
                                   #"Bilophila wadsworthia",
                                   #"Faecalibacterium prausnitzii"
                                   #"Parabacteroides distasonis",
                                   #"Paraprevotella xylaniphila",
                                   #"Coprobacter secundus",
                                   #"Coprobacter fastidiosus",
                                   #"Barnesiella intestinihominis"
                                   )) {                                
                        return(T)
                      } else {
                        return(F)
                      }
                    })) %>%
                    mutate(species = ifelse(!kickBool, "", species)),
                  aes(x = `PCo 1`, y = `PCo 2`, label = species), color = 'black', size = 3.5, min.segment.length = unit(0.5, 'lines'), max.overlaps = Inf) +
  theme_classic() +
  scale_color_manual(values = phylum_color_map) +
  theme(legend.text = element_text(size = 10)), filename = here('figures', "revisions", "Fig2_C_new.pdf"), width = 5.75, height = 3.5)


## Addition: Get motu-wise z-scores and scatter with PCs
library(readxl)
library(ggembl)
completed_substrate_annotations <- read_xlsx(here("data", "Glycan_Annotations", "20250219_Table_S1_incl_dbCAN3_annotations.xlsx"))
glycan_annotations_final_cleaned <- completed_substrate_annotations %>% select(c(Family:Subfamily,ORIGIN:FUNCTION_AT_DESTINATION_3, Glycan_annotation))
cazyAnnots <- glycan_annotations_final_cleaned

substrateNsBymOTU <- almeidaCAZy %>% left_join(cazyAnnots %>%
                                                        select(Subfamily,
                                                               FUNCTION_AT_DESTINATION_1) %>%
                                                        rename(cazy_family = Subfamily) %>%
                                                        separate_rows(FUNCTION_AT_DESTINATION_1, sep = ','), by = 'cazy_family') %>%
  group_by(genome, mOTU_ID, FUNCTION_AT_DESTINATION_1) %>%
  tally()

substrateNsBymOTU <- substrateNsBymOTU %>%
  group_by(mOTU_ID, FUNCTION_AT_DESTINATION_1) %>%
  summarize(nSD = sd(n),
            n = mean(n)) %>%
  filter(!is.na(FUNCTION_AT_DESTINATION_1)) %>%
  filter(!FUNCTION_AT_DESTINATION_1 == "Other") %>%
  filter(!is.na(mOTU_ID)) %>%
  replace(is.na(.), 0) %>% 
  pivot_wider(id_cols = mOTU_ID, names_from = FUNCTION_AT_DESTINATION_1, values_from = n, values_fill = 0)

# z-scores
substrateNsBymOTU[, !colnames(substrateNsBymOTU) == "mOTU_ID"] <- apply(substrateNsBymOTU[, !colnames(substrateNsBymOTU) == "mOTU_ID"], 2, scale)
substrateNsBymOTU <- substrateNsBymOTU %>%
  inner_join(pcoa %>% ungroup() %>% select(mOTU_ID)) %>%


pcoa_o <- pcoa %>%
  inner_join(substrateNsBymOTU, by = 'mOTU_ID')

plots <- list()
for (substrate in c(
  "DF",
  "Mucin",
  "GAG"
)) {
  tmp <- pcoa_o %>%
      mutate(substrate = get(substrate)) %>%
      mutate(`Substrate enrichment` = case_when(
        substrate > 3 ~ 3,
        substrate < -3 ~ -3,
        TRUE ~ substrate
      )) 
  p <- ggplot(tmp) +
    geom_point(aes(x = `PCo 1`, y = `PCo 2`, color = `Substrate enrichment`), alpha = 0.2) +
    theme_presentation() +
      scale_color_gradient2(
        low = "#0000FF",    # Blue for negative values
        mid = "#FFFFFF",    # White for midpoint (0)
        high = "#FF0000",   # Red for positive values
        midpoint = 0,       # Set the midpoint at 0
        limits = c(-3, 3)   # Cap the scale at -3 and 3
    ) +
    ggtitle(substrate)
    NULL
  plots[[length(plots) + 1]] <- p
}

ggsave(
  plot = wrap_plots(plots) +
    plot_layout(ncol = 3, guides = "collect"),
  filename = here('figures', "revisions", "Fig2_C_new_substrate_enrichments.pdf"), width = 8.5, height = 2.45
)


# cor_plots <- list()
# pcoa_o <- pcoa
# colnames(pcoa_o)[colnames(pcoa_o) %in% c('PCo 1', "PCo 2")] <- c("PCo_1", "PCo_2")
# pc_name_map <- c("PCo_1" = "PCo 1", "PCo_2" = "PCo 2")
# for (category in (substrateNsBymOTU %>% ungroup() %>% select(-c(mOTU_ID, Unknown)) %>% colnames())) {
#   for (pcname in c("PCo_1", "PCo_2")) {
#       tmp <- pcoa_o %>%
#         ungroup() %>%
#         select(mOTU_ID, {pcname}) %>%
#         left_join(substrateNsBymOTU %>% select(mOTU_ID, all_of(category)), by = 'mOTU_ID')
#       p <- ggplot() +
#         geom_point(data = tmp, aes_string(x = pcname, y = category), alpha = 0.2) +
#         theme_presentation()  +
#         xlab(pc_name_map[[pcname]]) +
#         ylab(category) +
#         theme(
#           axis.text.x = element_text(size = 6),
#           axis.text.y = element_text(size = 6),
#           axis.title.x = element_text(size = 7),
#           axis.title.y = element_text(size = 7),
#         ) +
#         NULL
#       # get cor
#       basecor <- stats::cor(tmp %>% select(all_of(c(pcname,category))) %>% as.data.frame() %>% as.matrix(), method = 'spearman')[1,2]
#       if (basecor > 0) {
#           p <- p + annotate(
#               'text',
#               x = -Inf, y = Inf,                     # Top-right corner
#               label = str_c("r = ", round(basecor, 2)),
#               size = 2.5,
#               hjust = -0.1, vjust = 1.1             # Adjust alignment to keep it inside the plot
#           )
#       } else {
#           p <- p + annotate(
#               'text',
#               x = Inf, y = Inf,                    # Top-left corner
#               label = str_c("r = ", round(basecor, 2)),
#               size = 2.5,
#               hjust = 1.1, vjust = 1.1            # Adjust alignment to keep it inside the plot
#           )
#       }
#       cor_plots[[length(cor_plots) + 1]] <- 
#       list(
#         plot = p,
#         cor = basecor,
#         category = category,
#         pcname = pcname
#       )
#   }
# }
# cor_plots <- enframe(cor_plots) %>%
#   mutate(
#     cor = map_dbl(value, \(x) x$cor),
#     category = map_chr(value, \(x) x$category),
#     pcname = map_chr(value, \(x) x$pcname),
#     plot = map(value, \(x) x$plot)
#   ) %>%
#   arrange(desc(abs(cor))) %>%
#   head(6) 

# library(patchwork)
# ggsave(
#   plot = wrap_plots(cor_plots$plot) +
#     plot_layout(ncol = 5, guides = "collect"),
#   filename = here('figures', "revisions", "Fig2_C_new_substrate_correlations.pdf"), width = 4, height = 2.75
# )


