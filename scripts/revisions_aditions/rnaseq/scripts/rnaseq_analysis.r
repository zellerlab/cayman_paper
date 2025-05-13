library(tidyverse)
library(patchwork)
library(tximport)
library(ggembl)
library(DESeq2)
library(readxl)


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
						tx2gene = tx,
						countsFromAbundance = "lengthScaledTPM")		
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

alldata_with_deseq2 <- alldata_with_deseq2 %>%
	left_join(cayman, by = 'strain') %>%
	mutate(deseq2_results = map2(deseq2_results, cayman_annotations, \(x, y) {
		if (!all(y$locus_tag %in% x$locus_tag)) {
			warning('fuck')
		}
		x <- left_join(x, y, by = 'locus_tag') %>%
			as_tibble()
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
	select(strain, baseMean, log2FoldChange, padj, gene_name, locus_tag, family, FUNCTION_AT_DESTINATION_1, mucin_related) %>%
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

tmp %>%
	mutate(strain = factor(strain, levels = unique(strain))) %>%
	filter(group == "sign. upregulated") %>%
	complete(family, strain, fill = list(log2FoldChange = NA)) %>%
	pivot_wider(
		id_cols = family, names_from = strain, values_from = log2FoldChange, values_fill = NA, values_fn = list(log2FoldChange = mean)
	)

#
