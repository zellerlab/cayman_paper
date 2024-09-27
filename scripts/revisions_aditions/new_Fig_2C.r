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
                                   "Hungatella hathewayi",
                                   "Dialister succinatiphilus",
                                   "Veillonella atypica",
                                   "Eggerthella lenta",
                                   "Bifidobacterium longum",
                                   "Bilophila wadsworthia",
                                   "Faecalibacterium prausnitzii"
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
                  aes(x = `PCo 1`, y = `PCo 2`, label = species, color = phylum), size = 3.5, min.segment.length = unit(0.5, 'lines'), max.overlaps = Inf) +
  theme_classic() +
  scale_color_manual(values = phylum_color_map) +
  theme(legend.text = element_text(size = 10)), filename = here('figures', "revisions", "Fig2_C_new.pdf"), width = 5.75, height = 3.5)
