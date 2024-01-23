library(here)

source(here('scripts', 'utils.r'))

# almeida_cazy_annotations_aggregated.tsv is generated and written at the beginning of fig_2.r
almeidaCAZyAggregated <- read_tsv(here('data', "almeida_cazy_annotations_aggregated.tsv"))
tree.filtered <- read.tree(here('data', 'tree.genus.ncbi.filtered.nwk'))

tmp <- almeidaCAZyAggregated %>%
  filter(genus %in% tree.filtered$tip.label) %>%
  select(-c(numberCAZyGenesSD, numberUniqueCAZyFamiliesSD, N)) %>%
  pivot_longer(-c(genus, familyGroup)) %>%
  rename(category = name)
k <- tmp %>% group_by(genus) %>% filter(category == 'numberCAZyGenes') %>% summarize(m = sum(value)) %>% arrange(desc(m)) %>% pull(genus) %>% head(50)
tmp <- tmp %>%
  filter(genus %in% k)

tmp$genus <- factor(tmp$genus, levels = k)
tmp$familyGroup <- factor(tmp$familyGroup, levels = rev(cazyFamilyOrder))

tmp <- tmp %>%
  rename(`copy number` = value)

ppp <- ggplot() +
  geom_bar(data = tmp, aes(x = genus, y = `copy number`, fill = familyGroup), stat = 'identity') +
  facet_grid(category ~ ., scales = 'free_y') +
  theme_classic() +
  scale_fill_manual(values = rev(familyColors)) +
  theme(axis.text.x = element_text(angle = 60, hjus = 1))

ggsave(plot = ppp, filename = here('figures', "Extended_Data_Fig1_genus_CAZy_richness.pdf"), width = 8, height = 5)