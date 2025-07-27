library(tidyverse)
library(patchwork)
library(tximport)
library(ggembl)
library(DESeq2)
library(readxl)
library(here)


# Total run time: ~5 minutes

rem_outlying_control_rep_in_Hhathewayi <- TRUE

# Load and consolidate salmon TMP data
p <- list.files('/g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/rnaseq/results/rnaseq_data_zusammengesucht_quantsf', full.names = T)
# Remove REP1 from control, H_hathewayi
if (rem_outlying_control_rep_in_Hhathewayi) {
	rem_bool <- str_detect(p, "CONTROL_REP1___H_hathewayi")
	p <- p[!rem_bool]
}

alldata <- map(p, \(x) {
	x <- read_tsv(x)
	return(x)
})

alldata <- enframe(alldata) %>%
	mutate(strain = p) %>%
	mutate(path_to_quantsf = strain) %>%
	mutate(strain = str_split_fixed(strain, "/", n = 20)[, 12]) %>%
	mutate(strain = str_replace(strain, "__salmon.merged.gene_tpm.tsv", "")) %>%
	#select(-name) %>%
	dplyr::rename(
		data = value) %>%
	mutate(controlrep = str_split_fixed(strain, "___", n = 3)[, 1]) %>%
	mutate(strain = str_split_fixed(strain, "___", n = 3)[, 2]) %>%
	arrange(desc(strain))

# Get TPM data and visualize in PCA

tpmdata <- alldata %>%
	group_by(strain) %>%
	nest() %>%
	mutate(data = map(data, \(x){
		#browser()
		x <- x %>%
			mutate(dataa = map2(data, controlrep, \(d, co) {
				
				return(d %>%
					mutate(controlrep = co) %>%
					select(Name, TPM, controlrep) %>%
					pivot_wider(id_cols = Name, names_from = controlrep, values_from = TPM))
			}))
		x <- purrr::reduce(x$dataa, full_join, by = "Name") %>%
			dplyr::rename(gene_name = Name)
		x <- x[rowSums(x[, -1]) > 0, ]		
		return(x)
	}))


# runs for 20 seconds
tpmdata <- tpmdata %>%
	mutate(pcoa_data = map(data, \(x) {
		x <- x %>%
			select(-gene_name) %>%
			t() %>%
			as.data.frame()
		x <- log10(x + 1)
		# pca, on log10-transformed euclidean distances, with feature scaling!
		pca_result <- prcomp(x, center = TRUE, scale. = TRUE)
		variance_explained <- summary(pca_result)$importance[2, 1:2]  # First two dimensions		
		return(list(pca_result, variance_explained))
	})) %>%
	mutate(variance_explained_pc1 = map_dbl(pcoa_data, \(x) x[[2]][1])) %>%
	mutate(variance_explained_pc2 = map_dbl(pcoa_data, \(x) x[[2]][2]))

tpmdata <- tpmdata %>%
	mutate(pcoa_plot = map2(pcoa_data, strain, \(x, s) {
		pca_result <- x[[1]]
		pca_data <- as.data.frame(pca_result$x)
		pca_data$sample <- rownames(pca_data)
		pca_data$treatment <- str_split_fixed(pca_data$sample, "_", 2)[, 1]
		pca_data$treatment <- case_when(
			pca_data$treatment == "TREATMENT" ~ "medium + mucin",
			pca_data$treatment == "CONTROL" ~ "medium")
		pca_data$replicate <- str_split_fixed(pca_data$sample, "_", 2)[, 2]
		pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = treatment, shape = replicate)) +
			geom_point() +
			xlab(paste0("PC1 (", round(x[[2]][1] * 100, 2), "%)")) +
			ylab(paste0("PC2 (", round(x[[2]][2] * 100, 2), "%)")) +
			theme_presentation() +
			ggtitle(s)
		return(pca_plot)
	}))

ggsave(plot = wrap_plots(tpmdata$pcoa_plot) + plot_layout(guides = 'collect'), filename = "../plots/pcoas_tpmlog10scaled.pdf", width = 12, height = 6)

# load tx2gene data
p_tx2gene <- list.files('/g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/rnaseq/results/rnaseq_data_zusammengesucht_tx2gene', full.names = T)
tx2gene <- map(p_tx2gene, \(x) {
	x <- read_tsv(x)
	return(x)
})
tx2gene <- enframe(tx2gene) %>%
	mutate(strain = p_tx2gene) %>%
	mutate(strain = str_split_fixed(strain, "/", n = 20)[, 12]) %>%
	mutate(strain = str_replace(strain, "___tx2gene.tsv", "")) %>%
	select(-name) %>%
	dplyr::rename(tx2gene = value)

# This whole block runs for around 2-3 minutes
alldata_with_deseq2 <- alldata %>%
	select(path_to_quantsf, controlrep, strain) %>%
	group_by(strain) %>%
	nest() %>% 
	left_join(tx2gene, by = 'strain') %>%
	mutate(deseq2_raw = map2(data, tx2gene, \(da, tx) {
		da <- da %>%
			mutate(control_or_treatment = str_split_fixed(controlrep, "_", 2)[, 1])
		da <- as.data.frame(da)
		na <- da$path_to_quantsf
		names(na) <- da$controlrep
		tx_object <- tximport(na, 
						type = "salmon", 
						tx2gene = tx
		)
		rownames(da) <- colnames(tx_object$counts)

		sample_info <- data.frame(condition = da$control_or_treatment)
		rownames(sample_info) <- colnames(tx_object$counts)

		dds <- DESeqDataSetFromTximport(tx_object, 
									colData = sample_info, 
									design = ~ condition)
		return(dds)
	})) %>%
	mutate(
		deseq2_results = map(deseq2_raw, \(x) {
			x <- DESeq(x)
			res <- results(x, contrast = c("condition", "TREATMENT", "CONTROL"))
			res <- as.data.frame(res)
			res$gene_name <- rownames(res)
			return(res)
		})
	)


# load gtf to be able to map to cayman
p_gtf <- list.files('/g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/rnaseq/results/rnaseq_data_zusammengesucht_gtf', full.names = T)
gtf <- map(p_gtf, \(x) {
	x <- read_tsv(x, skip = 3, comment = "#", col_names = F) %>%
		filter(X3 == "gene") %>%
		mutate(gene_name = str_split_fixed(X9, ";", 2)[, 1]) %>%
		mutate(gene_name = str_replace(gene_name, 'gene_id \"', '')) %>%
		mutate(gene_name = str_replace(gene_name, '\"', '')) %>%
		mutate(locus_tag = str_replace(X9, '.*locus_tag \"', '')) %>%
		mutate(locus_tag = str_replace(locus_tag, '\"', '')) %>%
		mutate(locus_tag = str_replace(locus_tag, ';.*', '')) %>%
		select(gene_name, locus_tag) %>% 
		distinct()
	return(x)
})
gtf <- enframe(gtf) %>%
	mutate(strain = p_gtf) %>%
	mutate(strain = str_split_fixed(strain, "/", n = 20)[, 12]) %>%
	mutate(strain = str_replace(strain, ".gtf", "")) %>%
	select(-name) %>%
	dplyr::rename(gtf = value)

alldata_with_deseq2 <- alldata_with_deseq2 %>%
	left_join(gtf, by = 'strain') %>%
	mutate(deseq2_results = map2(deseq2_results, gtf, \(x, y) {
		dim_before <- dim(x)[1]
		x <- x %>%
			inner_join(y, by = "gene_name")
		if(dim_before != dim(x)[1]) {
			warning(paste0("Number of genes in DESeq2 results changed from ", dim_before, " to ", dim(x)[1], " after joining with gtf."))
		}
		return(x)
	}))

# Load cazy annotations and mucin annotations and join.
p_cayman <- list.files('/g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/rnaseq/results/cayman_output', full.names = T)
cayman <- map(p_cayman, \(x) {
	x <- read_csv(x) %>%
		dplyr::rename(locus_tag = sequenceID)
	return(x)
})
cayman <- enframe(cayman) %>%
	mutate(strain = p_cayman) %>%
	mutate(strain = str_split_fixed(strain, "/", n = 20)[, 12]) %>%
	mutate(strain = str_replace(strain, ".cayman", "")) %>%
	select(-name) %>%
	dplyr::rename(cayman_annotations = value)

p_lociinfo <- list.files('/g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/rnaseq/data', full.names = T)
p_lociinfo <- p_lociinfo[str_detect(p_lociinfo, "loci.info")]
lociinfo <- map(p_lociinfo, \(x) {
	x <- read_tsv(x, col_names = F) %>%
		select(X1, X2, X5, X9, X12, X13) %>%
		mutate(X9_raw = X9) %>%
		filter(!str_detect(X2, "Dbxref")) %>%
		#filter(!is.na(X13)) %>%
		mutate(
			X1 = str_replace(X1, "gene_id \"", ""),
			X1 = str_replace(X1, "\"", ""),
			X2 = str_replace(X2, "transcript_id \"", ""),
			X2 = str_replace(X2, "\"", ""),
			X5 = str_replace(X5, "Name \"", ""),
			X5 = str_replace(X5, "\"", ""),
			X9 = str_replace(X9, "gene \"", ""),
			X9 = str_replace(X9, "\"", ""),
			X12 = str_replace(X12, "product \"", ""),
			X12 = str_replace(X12, "\"\t protein_id.*", ""),
			X13 = str_replace(X13, "product \"", ""),
			X13 = str_replace(X13, "\"\t protein_id.*", ""),			
		) %>%
		mutate(X12 = ifelse(str_detect(X12, 'locus_tag') | str_detect(X12, 'transl_table'), X13, X12)) %>%
		select(-X13) %>%
		mutate(X12 = str_replace(X12, "\"", ""))  %>%
		distinct()
	return(x)
})

lociinfo <- enframe(lociinfo) %>%
	mutate(strain = p_lociinfo) %>%
	mutate(strain = str_split_fixed(strain, "/", n = 20)[, 11]) %>%
	mutate(strain = str_replace(strain, "_from_agat.loci.info", "")) %>%
	select(-name) %>%
	dplyr::rename(lociinfoo = value)

alldata_with_deseq2 <- alldata_with_deseq2 %>%
	left_join(cayman, by = 'strain') %>%
	left_join(lociinfo, by = 'strain') %>%
	mutate(deseq2_results = pmap(list(deseq2_results, cayman_annotations, lociinfoo), \(x, y, z) {
		dim_y_before <- dim(y)[1]
		y <- inner_join(y, z, by = c("locus_tag" = "X5"))
		dim_y_after <- dim(y)[1]
		stopifnot(dim_y_before == dim_y_after)
		x <- left_join(x, y, by = c('gene_name' = "X1")) %>%
			as_tibble()
		stopifnot((x %>% filter(!is.na(family)) %>% dim() %>% magrittr::extract2(1)) == dim_y_after)
		return(x)		
	}))

# load annotations
annots <- read_xlsx('/g/scb/zeller/karcher/cayman_paper/data/Glycan_Annotations/20250219_Table_S1_incl_dbCAN3_annotations.xlsx') %>%
	select(Subfamily, FUNCTION_AT_DESTINATION_1)
alldata_with_deseq2_final <- alldata_with_deseq2 %>%
	mutate(deseq2_results = map(deseq2_results, \(x) {
		return(
			x %>% 
				left_join(annots, by = c("family" = "Subfamily")) %>% 
				mutate(mucin_related = case_when(
					is.na(FUNCTION_AT_DESTINATION_1) ~ "No CAZy annotation",
					str_detect(FUNCTION_AT_DESTINATION_1, "mucin") | str_detect(FUNCTION_AT_DESTINATION_1, "Mucin") ~ "Yes", 
					TRUE ~ "No")
				) %>%
	mutate(mucin_related = factor(mucin_related, levels = c("Yes", "No", "No CAZy annotation")))
				)
	})) %>%
	mutate(deseq2_results = map(deseq2_results, \(x) {
		return(x %>% arrange(padj))
	}))

# Get volcano plot
alldata_with_deseq2_final <- alldata_with_deseq2_final %>%
	mutate(VP = map2(deseq2_results, strain, \(x, ss) {
		p <- ggplot() +
			geom_point(data = x %>% filter(!mucin_related == "Yes"), aes(x = log2FoldChange, y = -log10(padj)), color = 'grey', alpha = 0.5) +
			geom_point(data = x %>% filter(mucin_related == "Yes"), aes(x = log2FoldChange, y = -log10(padj)), color = 'red', alpha = 0.75) +
			geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
			#geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
			theme_presentation() +
			xlab("log2 Fold Change") +
			ylab("-log10 adjusted p-value") +
			ggtitle(str_c(ss, "\nlog2FC > 0 -> expressed more highly in +mucin cond.", collapse = ""))
	}))

alldata_with_deseq2_final <- alldata_with_deseq2_final %>%
	mutate(VP_cut = map2(deseq2_results, strain, \(x, ss) {
		x <- x %>% filter(padj > 1E-50)
		p <- ggplot() +
			geom_point(data = x %>% filter(!mucin_related == "Yes"), aes(x = log2FoldChange, y = -log10(padj)), color = 'grey', alpha = 0.5) +
			geom_point(data = x %>% filter(mucin_related == "Yes"), aes(x = log2FoldChange, y = -log10(padj)), color = 'red', alpha = 0.75) +
			geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
			#geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
			theme_presentation() +
			xlab("log2 Fold Change") +
			ylab("-log10 adjusted p-value") +
			ggtitle(str_c(ss, "\nlog2FC > 0 -> expressed more highly in +mucin cond.\nOnly genes with padj > 1E-50", collapse = ""))
	}))


# save plots
devnull <- pmap(list(alldata_with_deseq2_final$VP, alldata_with_deseq2_final$VP_cut, alldata_with_deseq2_final$strain), \(x, x_cut, ss) {
				ggsave(plot = x + x_cut + plot_layout(ncol = 1), filename = paste0("../plots/volcano_", ss, ".pdf"), width = 6, height = 10)
				return(x)
			})

## Dig into the top overexpressed genes per strain by hand, for Georg (all hail the king)
tmp <- alldata_with_deseq2_final %>%
	ungroup() %>%
	select(strain, deseq2_results) %>%
	unnest(cols = c(deseq2_results)) %>%
	select(strain, baseMean, log2FoldChange, padj, gene_name, locus_tag.x, locus_tag.y, X9, X12, family, FUNCTION_AT_DESTINATION_1, mucin_related) %>%
	filter(mucin_related == "Yes") %>%
	mutate(group = case_when(
		log2FoldChange > 0 & padj < 0.05 ~ "sign. upregulated",
		log2FoldChange < 0 & padj < 0.05 ~ "sign. downregulated",
		TRUE ~ "not_significant"
	))

p <- tmp %>%
	group_by(strain, group) %>%
	tally() %>%
	mutate(group = factor(group, levels = c("sign. upregulated", "sign. downregulated", "not_significant"))) %>%
	ggplot() +
	geom_bar(aes(x = group, y = n, fill = group), stat = "identity") +
	facet_wrap(~ strain, ncol = 2) +
	theme_presentation() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(
	plot = p,
	filename = "../plots/overexpressed_genes_per_strain.pdf",
	width = 7.5,
	height = 7.5
)


families_2B <- c("GH2", "GH92", "GH20", "GH31", "GH29", "GH97", "GH95", "CBM32", "GH36", "GH35", 
  "GH33", "GH42", "GH130", "GH18", "GH16", "GH109", "CBM51", "GH110", "GH84", 
  "GH27", "GH89", "GH123", "GH85", "GH136", "GH112")

dat <- alldata_with_deseq2_final %>%
	ungroup() %>%
	select(strain, deseq2_results) %>%
	unnest(cols = c(deseq2_results)) %>%
	filter(baseMean != 0) %>%
	select(strain, baseMean, log2FoldChange, padj, gene_name, locus_tag.x, locus_tag.y, X9, X12, family, FUNCTION_AT_DESTINATION_1, mucin_related) %>%
	#filter(padj < 0.05) %>%
	inner_join(data.frame(family = families_2B)) %>%
	mutate(family = factor(family, levels = families_2B)) %>%
	mutate(
		`Expression\ndifferences` = case_when(
			log2FoldChange > 0 & padj < 0.05 ~ "upregulated",
			log2FoldChange < 0 & padj < 0.05 ~ "downregulated",
			TRUE ~ "not diff. regulated"
		)
	) %>%
	mutate(
		`Expression\ndifferences` = factor(`Expression\ndifferences`, levels = c("upregulated", "downregulated", "not diff. regulated"))
	) %>%
	mutate(
		strain = case_when(
			strain == "H_hathewayi_DSM13479" ~ "Hungatella hathewayi (DSM13479)",
			strain == "E_tayi_DSM26961" ~ "Eisenbergiella tayi (DSM26961)"
		)
	)

t_test_res <- dat %>%
	group_by(strain) %>%
	summarise(
		t_test_p = if (n() > 2) t.test(log2FoldChange, alternative = "greater")$p.value else NA_real_
	) %>%
	ungroup() %>%
	mutate(t_test_p_adj = p.adjust(t_test_p, method = "BH"))

library(ggrepel)
p <- ggplot() + 
	geom_hline(
		yintercept = 0, linetype = "solid", color = "black", alpha = 0.2
		) +
	#geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0, dodge.width = 0.5)) +
	geom_point(data = dat, aes(x = family, y = log2FoldChange, alpha = `Expression\ndifferences`)) +
	geom_point(data = dat %>% anti_join(tibble(`Expression\ndifferences` = "not diff. regulated")), aes(x = family, y = log2FoldChange, color = `Expression\ndifferences`), shape = 1, size = 3) +
	geom_text_repel(data = dat %>% anti_join(tibble(`Expression\ndifferences` = "not diff. regulated")), aes(x = family, y = log2FoldChange, label = locus_tag.x), size = 1.5, segment.size = 0.25) +
	theme_presentation() + 
	theme(
		axis.text.x = element_text(angle = 45,  hjust = 1, size = 6),
		axis.text.y = element_text(size = 7)
		) +
	scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "not diff. regulated" = "grey")) +
	scale_alpha_manual(values = c("upregulated" = 0.175, "downregulated" = 0.175, "not diff. regulated" = 0.175)) +
	facet_wrap(.~strain, ncol = 1) +
	guides(alpha = "none") +
	ylab("Fold change (log2)") + 
    theme(
        legend.position = "bottom",          # Place legend at the bottom
        legend.direction = "horizontal",     # Make legend horizontal
        legend.box = "vertical"              # Stack legend entries vertically
    ) +
    guides(color = guide_legend(ncol = 1)) +
	scale_x_discrete(drop = FALSE) +
	ylim(c(-5.5, 8.5))

ggsave(
	plot = p,
	filename = here("scripts", "revisions_aditions", "rnaseq", "plots", "FCs_per_gene_and_family.pdf"),
	width = 3.25,
	height = 3.25
)	

# Also try random effects modelling
library(lme4)
library(lmerTest)
llm_dat <- alldata_with_deseq2_final %>%
	ungroup() %>%
	select(strain, deseq2_results) %>%
	unnest(cols = c(deseq2_results)) %>%
	filter(baseMean != 0) %>%
	select(strain, family, X2) %>%
	#filter(padj < 0.05) %>%
	inner_join(data.frame(family = families_2B)) %>%
	mutate(family = factor(family, levels = families_2B)) %>%
	#mutate(gene_name = str_replace(gene_name, "-gene", '-rna'))  %>%
	filter(!is.na(family)) %>%
	left_join(
		tpmdata %>% 
			select(strain, data) %>% 
			unnest()
		, by = c("strain", 'X2' = 'gene_name')	
	) %>%
	pivot_longer(-c(strain, X2, family), names_to = "rep", values_to = "TPM") %>%
	dplyr::rename(gene_name = X2) %>%
	# Some notes
	# For H. hathewayi, control_rep3 is all NAs (we kicked it out). 
	# TODO: Filter this once data is in long (further down)
	filter(!is.na(TPM)) %>%
	mutate(sample_type = ifelse(str_detect(rep, "TREATMENT"), 'With mucin', 'Without mucin')) %>%
	group_by(strain, family, gene_name) %>%
	mutate(TPM_scaled = scale(TPM)[, 1]) 

llm_res <- llm_dat %>%
	group_by(strain) %>%
	nest() %>%
	mutate(lme_res_version1 = map(
		data, \(x) {
			mm <- x$TPM[x$TPM > 0]
			mm <- min(mm)

			x <- x %>%
				lmer(scale(log10(TPM + mm)) ~ sample_type + (1 | gene_name), data = .)
				#lmer(TPM ~ sample_type + (1 | gene_name), data = .)
			return(x)
		}
	)) %>%
	mutate(lme_res_version2 = map(
		data, \(x) {
			mm <- x$TPM[x$TPM > 0]
			mm <- min(mm)	

			x <- x %>%
				lmer(scale(log10(TPM + mm)) ~ sample_type + (1 | family) + (1 | gene_name), data = .)
				#lmer(TPM ~ sample_type + (1 | gene_name), data = .)
			return(x)
		}
	))

for (strain in llm_res$strain) {
	for (version in c("version1", "version2")) {
		print(str_c("strain: ", strain, " (", version, ")"))
		print(summary(llm_res[[str_c("lme_res_", version)]][[which(llm_res$strain == strain)]])$coefficients)
	}
}


# expression ~ mucin_condition + (1 | CAzyme_family) + (1 | gene)


