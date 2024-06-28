suppressPackageStartupMessages({library(ggplot2)
        library(dplyr)
        library(scales)
        library(tidyr)
    	library(tibble)
        library(egg)
		library(data.table)})

main <- function() {
	## process data 
	res_file = snakemake@input$variantOverlap
	p_threshold = snakemake@params$p_threshold %>% as.numeric()
	out_enr = snakemake@output$outFile_enrichment
	out_recall = snakemake@output$outFile_recall

	res = fread(res_file, sep="\t", header=TRUE)
	res = dplyr::filter(res, is.finite(enrichment), is.finite(recall), nVariantsTotal>20)
	res$enhMb = res$bpEnhancers/1e6

	# cluster (both based on enrichment for consistency)
	M = dplyr::select(res, biosample, trait, enrichment) %>%
		pivot_wider(names_from=trait, values_from=enrichment) %>% column_to_rownames("biosample")
	M[is.na(M)] = 0
	
	trait_dist = dist(1-cor(M))
	trait_dist[is.na(trait_dist)] = 0
	order_traits = hclust(trait_dist, method = "ward.D2")$order
	traits_ordered = colnames(M)[order_traits]

	biosample_dist = dist(1-cor(t(M)))
	biosample_dist[!is.na(biosample_dist)] = 0
	order_biosamples = hclust(biosample_dist, method="ward.D2")$order
	biosamples_ordered = rownames(M)[order_biosamples]

	res$trait = factor(res$trait, levels=traits_ordered, ordered=TRUE)
	res$biosample = factor(res$biosample, levels=biosamples_ordered, ordered=TRUE)

	# plotting params
	enr_colors =  c("#f6eff7","#bdc9e1", "#67a9cf","#1c9099", "#016c59")
	recall_colors = c("#edf8fb","#b3cde3", "#8c96c6", "#8856a7", "#810f7c")
	bar_colors = c("#435369", "#96a0b3")
	na_color = "#ffffff"
	enr_lims = c(0, max(res$enrichment))
	recall_limits = c(0, max(res$recall))
	ht = ifelse(length(rownames(M))>50, 16, 8) # dependent on number of biosamples

	## plots
	# enrichment heatmap
	enr_sig = dplyr::filter(res, p_adjust_enr<p_threshold)
	enr_grid = ggplot(res, aes(x=trait, y=biosample, fill=enrichment)) +
		geom_tile() +
		geom_point(data=enr_sig, shape=8, size=0.25) + 
		scale_fill_gradientn(colors=enr_colors, oob=scales::squish, na.value=na_color, limits=enr_lims, name="Enrichment") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_blank(), axis.text.x = element_blank(), panel.border = element_blank(),
			legend.position="top", legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7))

	# recall heatmap
	rec_grid = ggplot(res, aes(x=trait, y=biosample, fill=recall)) +
		geom_tile() +
		scale_fill_gradientn(colors=recall_colors, oob=scales::squish, na.value=na_color, limits=recall_limits, name="Recall") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_blank(), axis.text.x = element_blank(),
			legend.position='top',  legend.direction='horizontal', legend.text=element_text(size=7), legend.title=element_text(size=7))

	# variants per trait
	n_var = dplyr::select(res, trait, nVariantsTotal) %>% distinct()
	var_count = ggplot(n_var, aes(x=trait, y=nVariantsTotal)) +
		geom_bar(stat="identity", width=0.5, fill=bar_colors[1]) +
		ylab("# variants per trait") + xlab("") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.text.x = element_text(angle=60, hjust=1))

	# enhancer size
	enh_size = dplyr::select(res, biosample, group, enhMb) %>% distinct()
	enh_mb = ggplot(enh_size, aes(x=biosample, y=enhMb, fill=group)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("Enhancer set size\n(Mb)") + xlab("") +
		scale_fill_manual(values=bar_colors) + 
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), axis.text.y = element_blank(), legend.position="None") + 
		coord_flip()
	
	# assemble=
	blank = ggplot() + theme_void()

	enr_assembled = egg::ggarrange(enr_grid, enh_mb, var_count, blank, nrow=2, ncol=2, heights=c(2, 0.2), widths=c(2, 0.3))
	rec_assembled = egg::ggarrange(rec_grid, enh_mb, var_count, blank, nrow=2, ncol=2, heights=c(2, 0.2), widths=c(2, 0.3))

	ggsave(out_enr, enr_assembled, width=10, height=ht)
	ggsave(out_recall, rec_assembled, width=10, height=ht)
}

main()
