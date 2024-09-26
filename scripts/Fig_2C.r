library(here)
library(ggrepel)
source(here('scripts', 'utils.r'))

# motus3.0_taxonomy <- read_tsv(here('data', "motus3.0_taxonomy.tsv"))
almeidaCAZy <- read_tsv(here('data', 'almeida_cazy_annotations.tsv'))
motus3.0_taxonomy <- read_tsv(here('data', "motus3.0_taxonomy.tsv"))
tree.filtered <- read.tree(here('data', "tree.genus.ncbi.filtered.nwk"))

uniqueFamiliesRaw <- almeidaCAZy %>%
  # We're only interested in Genera that we can name...
  filter(!str_detect(Genus, "gen[.]")) %>%
  filter(!str_detect(Genus, "incertae")) %>%
  group_by(Genus) %>%
  nest() %>%
#   mutate(Genus = map_chr(Genus, function(x) str_c(str_split(x, " ")[[1]][2:length(str_split(x, " ")[[1]])],
#                                                   sep = " ", collapse = " ")))
  mutate(Genus = str_split_fixed(Genus, " ", n=2)[, 2])
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

cazyFamilyPairwiseJaccards <- get_family_sharing_rate(uniqueFamiliesRaw, tree.filtered$tip.label, lev = "Genus")
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
                            column_to_rownames('taxon_1') %>%
                            as.matrix()))), k = 2)
colnames(pcoa) <- c("PCo 1", "PCo 2")
ggsave(plot = ggplot() + geom_point(data = pcoa %>%
                                     as.data.frame() %>%
                                     mutate(genus = rownames(.),
                                            #yOffset = `PCo 2` + 0.015) %>%
                                            yOffset = `PCo 2`) %>%           
                                     left_join(motus3.0_taxonomy %>%
                                                 mutate(Genus = str_split_fixed(Genus, " ", n = 2)[, 2]) %>%
                                                 mutate(Phylum = str_split_fixed(Phylum, " ", n = 2)[, 2]) %>%
                                                 select(Genus, Phylum) %>% distinct() %>%
                                                 rename(genus = Genus,
                                                        phylum = Phylum), by = 'genus') %>%
                                     # Collinsella has an ill-defined phylum. Fix.
                                     mutate(phylum = ifelse(genus == "Collinsella", "Actinobacteria", phylum)) %>%
                                     mutate(`PCo 2` = -1 * `PCo 2`) %>%
                                     mutate(`PCo 1` = -1 * `PCo 1`), aes(x = `PCo 1`, y = `PCo 2`, color = phylum)) +
  #geom_text(aes(label = genus, y = yOffset)) +
  geom_text_repel(data = pcoa %>%
                    as.data.frame() %>%
                    mutate(genus = rownames(.),
                           #yOffset = `PCo 2` + 0.015) %>%
                           yOffset = `PCo 2`) %>%           
                    left_join(motus3.0_taxonomy %>%
                                mutate(Genus = str_split_fixed(Genus, " ", n = 2)[, 2]) %>%
                                mutate(Phylum = str_split_fixed(Phylum, " ", n = 2)[, 2]) %>%
                                select(Genus, Phylum) %>% distinct() %>%
                                rename(genus = Genus,
                                       phylum = Phylum), by = 'genus') %>%
                    # Collinsella has an ill-defined phylum. Fix.
                    mutate(phylum = ifelse(genus == "Collinsella", "Actinobacteria", phylum)) %>%
                    mutate(kickBool = pmap_lgl(list(`PCo 1`, `PCo 2`, genus), function(po1, po2, g) {
                      # These are equivalent to the top 8 mucin-targetting taxa from 2B
                      if (g %in% c("Akkermansia",
                                   "Bacteroides",
                                   "Eisenbergiella",
                                   "Hungatella",
                                   "Parabacteroides",
                                   "Paraprevotella",
                                   "Coprobacter",
                                   "Barnesiella")) {                                
                        return(T)
                      } else {
                        return(F)
                      }
                    })) %>%
                    mutate(genus = ifelse(!kickBool, "", genus)) %>%
                    mutate(`PCo 2` = -1 * `PCo 2`) %>%
                    mutate(`PCo 1` = -1 * `PCo 1`),
                  aes(x = `PCo 1`, y = `PCo 2`, label = genus, color = phylum), size = 3.5, min.segment.length = unit(0, 'lines')) +
  theme_classic() +
  theme(legend.text = element_text(size = 10)), filename = here('figures', "Fig2_C.pdf"), width = 5.75, height = 3.5)