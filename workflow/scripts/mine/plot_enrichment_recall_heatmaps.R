suppressPackageStartupMessages({library(ggplot2)
        library(dplyr)
        library(scales)
        library(tidyr)
    	library(tibble)
        library(egg)
		library(data.table)})

cluster_traits_biosamples <- function(res){
	M = dplyr::select(res, biosample, trait, metric) %>%
	pivot_wider(names_from=trait, values_from=metric) %>% column_to_rownames("biosample")
	M[is.na(M)] = 0
	
	trait_dist = dist(1-cor(M))
	trait_dist[is.na(trait_dist)] = 0
	order_traits = hclust(trait_dist, method = "ward.D2")$order
	traits_ordered = colnames(M)[order_traits]

	biosample_dist = dist(1-cor(t(M)))
	biosample_dist[is.na(biosample_dist)] = 0
	order_biosamples = hclust(biosample_dist, method="ward.D2")$order
	biosamples_ordered = rownames(M)[order_biosamples]

	res$trait = factor(res$trait, levels=traits_ordered, ordered=TRUE)
	res$biosample = factor(res$biosample, levels=biosamples_ordered, ordered=TRUE)
	return(res)
}

plot_heatmap <- function(res, metric, metric_name, colors, max_value, out_file, p_threshold){
	# subset data
	res$metric = res[[metric]]
	
	# cluster
	res = cluster_traits_biosamples(res)

	# plotting params
	na_color = "#ffffff"
	text_color = "#ffffff"
	bar_colors = c("#435369", "#96a0b3")
	metric_lims = c(0, max_value)

	ht = ifelse(length(unique(res$biosample))>60, 12, 8) # dependent on number of biosamples
	
	# heatmap
	hm = ggplot(res, aes(x=trait, y=biosample, fill=metric)) +
		geom_tile() +
		scale_fill_gradientn(colors=colors, oob=scales::squish, na.value=na_color, limits=metric_lims, name=metric_name) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_blank(), axis.text.x = element_blank(), panel.border = element_blank(),
			legend.position="top", legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7))

	if (metric=="enrichment"){
		enr_sig = dplyr::filter(res, p_adjust_enr<p_threshold)
		hm = hm + geom_point(data=enr_sig, shape=8, size=0.25, color=text_color)
	}

	## margin plots
	# variants per trait
	n_var = dplyr::select(res, trait, traitGroup, nVariantsTotal) %>% distinct()
	var_max = 1000
	n_var$nVariantsTotal = pmin(n_var$nVariantsTotal, var_max)
	var_count = ggplot(n_var, aes(x=trait, y=(nVariantsTotal), fill=traitGroup)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("Variants / trait") + xlab("") +
		scale_fill_manual(values=bar_colors) + 
		theme_classic() + theme(axis.text = element_text(size = 7), axis.text.x = element_text(angle=60, hjust=1), legend.position="None")

	# enhancer size
	enh_size = dplyr::select(res, biosample, biosampleGroup, enhMb) %>% distinct()
	enh_mb = ggplot(enh_size, aes(x=biosample, y=(enhMb), fill=biosampleGroup)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("Enhancer\nset size) (Mb)") + xlab("") +
		scale_fill_manual(values=bar_colors) + 
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.text.y = element_blank(), legend.position="None") + 
		coord_flip()

	# combine with margin plots
	blank = ggplot() + theme_void()
	assembled = egg::ggarrange(hm, enh_mb, var_count, blank, nrow=2, ncol=2, heights=c(2, 0.2), widths=c(2, 0.3))
	ggsave(out_file, assembled, width=10, height=ht)
}


main <- function() {
	## process data 
	res_files <- c("/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results/2024_0911_hd_blood/scE2G_multiome/variant_overlap/enrichmentRecall.thresholded.traitByBiosample.tsv.gz",
		"/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results/2024_0829_pancreas/scE2G_multiome/variant_overlap/enrichmentRecall.thresholded.traitByBiosample.tsv.gz")
	max_enr <- 15
	max_recall <- 0.30
	p_threshold <- 0.05
	out_dir <- c("/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0823_GWAS_variant_interpretation")

	# plotting params
	enr_ovl_colors =  c("#f1eef6","#bdc9e1", "#74a9cf","#2b8cbe", "#045a8d") # PuBu
	recall_ovl_colors = c("#edf8fb","#b3cde3", "#8c96c6", "#8856a7", "#810f7c") # BuPu
	precision_colors = c("#edf8fb","#b2e2e2", "#66c2a4", "#2ca25f", "#006d2c")  # BuGn
	recall_link_colors = c("#ffffcc","#c2e699", "#78c679", "#31a354", "#006837") # YlGn

	res <- lapply(res_files, fread) %>% rbindlist()
	res <- dplyr::filter(res, is.finite(enrichment), is.finite(recall), !is.na(enrichment), !is.na(recall), nVariantsTotal>20,
		traitGroup == FALSE, biosampleGroup == FALSE)
	res$enhMb = res$bpEnhancers/1e6
	res$enhKb = res$bpEnhancers/1e3

	## plot enrichment heatmap
	plot_heatmap(res, "enrichment", "Enrichment", enr_ovl_colors, max_enr, file.path(out_dir, "GWAS_enrichment_heatmap_all.pdf"), p_threshold)

	## plot recall heatmap
	plot_heatmap(res, "recall", "Recall", recall_ovl_colors, max_recall,  file.path(out_dir, "GWAS_recall_heatmap_all.pdf"), p_threshold)
}

main()
