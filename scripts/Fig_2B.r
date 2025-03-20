library(here)
library(readxl)
source(here('scripts', 'utils.r'))

motus3.0_taxonomy <- read_tsv(here('data', "motus3.0_taxonomy.tsv"))
almeidaCAZy <- read_tsv(here('data', 'almeida_cazy_annotations.tsv'))
#cazyAnnots <- read_tsv(here('data', "20230609_glycan_annotations_cleaned_Nic.tsv"))
completed_substrate_annotations <- read_xlsx(here("data", "Glycan_Annotations/", "20250219_Table_S1_incl_dbCAN3_annotations.xlsx"))
glycan_annotations_final_cleaned <- completed_substrate_annotations %>% select(c(Family:Subfamily,ORIGIN:FUNCTION_AT_DESTINATION_3, Glycan_annotation))
cazyAnnots <- glycan_annotations_final_cleaned
tree.filtered <- read.tree(here('data', "tree.genus.ncbi.filtered.nwk"))

# This runs for a few minutes
mucin_degradation_pathway_abundances <- almeidaCAZy %>% 
  group_by(sequenceID, 
           genome, 
           cazy_family, 
           mOTU_ID, 
           Species,
           Genus) %>% 
  tally()  %>%  
  left_join(cazyAnnots, by = c('cazy_family' = 'Subfamily')) %>% 
  separate_rows(FUNCTION_AT_DESTINATION_1, sep = ",")

mucin_degradation_pathway_abundances <- mucin_degradation_pathway_abundances %>%
  filter(FUNCTION_AT_DESTINATION_1 == "Mucin") %>% 
  group_by(genome, 
           cazy_family) %>% 
  tally() %>%
  pivot_wider(id_cols = genome, 
              names_from = cazy_family, 
              values_from = n, 
              values_fill = 0) %>% 
  pivot_longer(-genome) %>% 
  rename(cazy_family = name, n = value) %>%
  left_join(almeidaCAZy %>%
              ungroup() %>% 
              select(genome, Genus, Species) %>%
              distinct(), by = "genome") %>%
  mutate(Genus = str_split_fixed(Genus, " ", n = 2)[, 2])

mucin_degradation_pathway_abundances <- mucin_degradation_pathway_abundances %>% 
  pivot_wider(id_cols = c(genome, Genus, Species), names_from = cazy_family, values_from = n, values_fill = 0) %>%
  group_by(Species) %>%
  pivot_longer(-c(genome, Genus, Species)) %>%
  rename(cazy_family = name, n = value) %>%
  group_by(Genus, cazy_family) %>%
  summarize(m = mean(n),
            e = sd(n),
            numberGenomes = length(Species)) %>%
  filter(!is.na(Genus)) %>%
  mutate(e = ifelse(is.na(e), 0 , e))

mucin_degradation_pathway_abundances <- mucin_degradation_pathway_abundances %>%
  filter(Genus %in% tree.filtered$tip.label)

mucin_degradation_pathway_abundances <- mucin_degradation_pathway_abundances %>%
  inner_join(mucin_degradation_pathway_abundances %>% group_by(Genus) %>% summarize(m = sum(m)) %>%
  arrange(desc(m)) %>%
  head(8) %>%
  select(Genus)) %>%
  left_join(motus3.0_taxonomy %>% 
  select(Phylum, Genus) %>% 
  distinct() %>% 
  mutate(Phylum = str_split_fixed(Phylum, " ", n = 2)[, 2]) %>% 
  mutate(Genus = str_split_fixed(Genus, " ", n = 2)[, 2]) %>% 
  distinct()) %>%
    mutate(Genus = map2_chr(Genus, numberGenomes, function(a, b) return(str_c(a, " (", b, ")", sep = "", collapse = "")))) %>%
    group_by(cazy_family) %>%
    filter(any(m > 1))

mucin_degradation_pathway_abundances %>%
  select(Genus, Phylum) %>%
  distinct() %>%
  arrange(Phylum)

l <- mucin_degradation_pathway_abundances %>%
  group_by(cazy_family) %>%
    summarize(m = mean(m)) %>%
    arrange(desc(m)) %>%
    filter(m != 0) %>%
  pull(cazy_family)
u <- length(unique(mucin_degradation_pathway_abundances$cazy_family))
w <- 0.8

mucin_degradation_pathway_abundances$Genus <- factor(mucin_degradation_pathway_abundances$Genus, levels = names(mucin_pathway_colors))

p <- ggplot(mucin_degradation_pathway_abundances %>% 
                   mutate(cazy_family = factor(cazy_family, levels = l)), aes(x = cazy_family, y = m, color = Genus, shape = Genus)) +
  #geom_jitter(size = 1.5, alpha = 0.5, height = 0, width = 0.1) +
  #geom_errorbar(aes(ymin = m - e, ymax = m + e), position='dodge', width = w, color = 'black', alpha = 0.5) +
  geom_linerange(aes(ymin = m - e, ymax = m + e), position=position_dodge(width = w), color = 'black', alpha = 0.25) +
  geom_point(size = 1.5, position = position_dodge(width = w)) +
  #geom_line(aes(group = Genus)) + 
  theme_classic() +
  scale_shape_manual(values = c(20,20,20,20,20,20,20,20)) +
  theme(axis.text.x = element_text(angle = 45, hjust =1)) +
  geom_vline(xintercept = (1:u)+0.5, alpha = 0.2, linetype = 'dashed') + 
  ylab("copy number (mean)") +
  xlab("CAZy family") + 
  #scale_y_log10() + 
  guides(color=guide_legend(nrow=4, byrow = TRUE),
         shape = guide_legend(nrow=4, byrow = TRUE)) +
  #scale_y_continuous(limits = c(0,10), breaks = c(0, 2, 4, 6, 8, 10)) +
  #scale_color_manual(values = rainbow(13)[c(1,5,8,9,11,12,13,6)]) +
  scale_color_manual(values = mucin_pathway_colors) +
  theme(plot.margin = margin(0.5, 0, 0, 0, "pt"),
        legend.position = 'bottom',
        legend.text = element_text(size = 10))

ggsave(plot = p, filename = here("figures", "Fig2_B.pdf"), width = 5.5, height = 5.5)