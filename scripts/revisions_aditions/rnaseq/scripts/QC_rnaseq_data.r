library(tidyverse)
library(ggembl)
library(patchwork)

# for f in $(find ../results/nfcore_rnaseq/ | grep star_summary); do ff=$(echo $f | rev | cut -d "/" -f 5 | rev); cat ${f} | sed "s%^%${ff}\t%"; done  > ../results/star_summary_all.tsv; 
# Then, manually remove x[1,1] field to make it readable

star_summary_all <- read_tsv("../results/star_summary_all.tsv", col_names = TRUE) %>%
	select(
		Tax, Sample, `Total reads`, `Aligned`, `Mismatch rate`, `Uniq aligned`
	) %>%
	rename(
		tax_raw = Tax,
		sample = Sample,
		total_reads_million = `Total reads`,
		aligned_percentage = `Aligned`,
		aligned_uniquely_percentage = `Uniq aligned`,
		mismatch_percentage = `Mismatch rate` # I know this says rate but looking at the multiQC plot, I see a percentage.
	) %>%
	filter(sample != "Sample") %>%	
	mutate(across(c(
		'aligned_percentage',
		'total_reads_million',
		'mismatch_percentage',
		"aligned_uniquely_percentage"
		), as.numeric)) %>%
	mutate(
		species = str_split_fixed(tax_raw, "_", 3)[, 1:2] %>% apply(1, paste0, collapse = "_"),
		strain = str_split_fixed(tax_raw, "_", 3)[, 3],
		replicate = str_split_fixed(sample, "_", 2)[, 2],
		treatment = str_split_fixed(sample, "_", 2)[, 1],
		total_reads = as.numeric(str_replace(total_reads_million, "M", "")) * 1e6,
	)  %>%
	mutate(
		treatment = case_when(
			treatment == "CONTROL" ~ "Without mucin",
			treatment == "TREATMENT" ~ "With mucin"
		)) %>%
	mutate(treatmant = factor(treatment, levels = c("Without mucin", "With mucin"))) %>%
	mutate(
		strain = case_when(
			species == "E_tayi" ~ "Eisenbergiella tayi (DMS26961)",
			species == "H_hathewayi" ~ "Hungatella hathewayi (DMS13479)",
		)
	) %>%
	mutate(
		species = case_when(
			species == "E_tayi" ~ "Eisenbergiella tayi",
			species == "H_hathewayi" ~ "Hungatella hathewayi",
		)
	)	

p0.9 <- ggplot() +
	geom_point(data = star_summary_all, aes(x = treatment, y = total_reads, shape = replicate)) +
	facet_grid(.~strain) +
	theme_presentation() +
	ylab("Total num reads") +
	theme(
		# Remove everything from the x axis
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.title.x = element_blank()

	) + 
	NULL


p1 <- ggplot() +
	geom_point(data = star_summary_all, aes(x = treatment, y = aligned_percentage, shape = replicate)) +
	facet_grid(.~strain) +
	theme_presentation() +
	ylab("% Aligned reads") +
	theme(
		# Remove everything from the x axis
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.title.x = element_blank()

	) + 
	NULL

p1.1 <- ggplot() +
	geom_point(data = star_summary_all, aes(x = treatment, y = aligned_uniquely_percentage, shape = replicate)) +
	facet_grid(.~strain) +
	theme_presentation() +
	ylab("% Aligned reads\n(uniquely)") +
	theme(
		# Remove everything from the x axis
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.title.x = element_blank()
	) + 
	NULL

p2 <- ggplot() +
	geom_point(data = star_summary_all, aes(x = treatment, y = mismatch_percentage, shape = replicate)) +
	facet_grid(.~strain) +
	theme_presentation() +
	ylab("% mismatches") +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1)
	) +
	#scale_y_log10() +
	# scale_y_log10() and force a tick at 0.1
	scale_y_log10(breaks = c(0.1, 0.5, 1, 1.5), limits = c(0.05, 2)) +
	NULL

ggsave(
	filename = "../plots/star_summary_all.pdf",
	plot = p0.9 + p1 + p1.1 + p2 + plot_layout(ncol = 1, guides = 'collect'),
	width = 6.75,
	height = 8
)
