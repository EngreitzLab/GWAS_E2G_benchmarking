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
	exc_all = dplyr::filter(n_var, trait!="ALL") 
	var_max = max(exc_all$nVariantsTotal)
	n_var$nVariantsTotal = pmin(n_var$nVariantsTotal, var_max)
	var_count = ggplot(n_var, aes(x=trait, y=log10(nVariantsTotal), fill=traitGroup)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("log10(variants per trait)") + xlab("") +
		scale_fill_manual(values=bar_colors) + 
		theme_classic() + theme(axis.text = element_text(size = 7), axis.text.x = element_text(angle=60, hjust=1), legend.position="None")

	# enhancer size
	enh_size = dplyr::select(res, biosample, biosampleGroup, enhKb) %>% distinct()
	enh_max = max(enh_size$enhKb)
	enh_size$enhKb = pmin(enh_size$enhKb, enh_max)
	enh_mb = ggplot(enh_size, aes(x=biosample, y=log10(enhKb), fill=biosampleGroup)) +
		geom_bar(stat="identity", width=0.5) +
		ylab("log10(enhancer\nset size) (kb)") + xlab("") +
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
	res_file = snakemake@input$variantOverlap
	ranges_file = snakemake@input$ranges
	fixed_scale = snakemake@params$fixedScale
	p_threshold = snakemake@params$p_threshold %>% as.numeric()
	out_enr = snakemake@output$outFile_enrichment
	out_recall = snakemake@output$outFile_recall

	# plotting params
	enr_ovl_colors =  c("#f1eef6","#bdc9e1", "#74a9cf","#2b8cbe", "#045a8d") # PuBu
	recall_ovl_colors = c("#edf8fb","#b3cde3", "#8c96c6", "#8856a7", "#810f7c") # BuPu
	precision_colors = c("#edf8fb","#b2e2e2", "#66c2a4", "#2ca25f", "#006d2c")  # BuGn
	recall_link_colors = c("#ffffcc","#c2e699", "#78c679", "#31a354", "#006837") # YlGn

	res = fread(res_file, sep="\t", header=TRUE)
	res = dplyr::filter(res, is.finite(enrichment), is.finite(recall), !is.na(enrichment), !is.na(recall), nVariantsTotal>20)
	res$enhMb = res$bpEnhancers/1e6
	res$enhKb = res$bpEnhancers/1e3

	# define max for color scales
	if (fixed_scale){
		ranges = fread(ranges_file, sep="\t")

		# get max enr
		enr = dplyr::filter(ranges, metric=="enrichment_overlap", min_max == "max")
		max_enr = quantile(enr$value, c(0.75))[1] # use 75th percentile max 

		# get max recall
		recall = dplyr::filter(ranges, metric=="recall_overlap", min_max == "max")
		max_recall = max(recall$value)

	} else {
		res_temp = dplyr::filter(res, nVariantsTotal>25, nCommonVariantsOverlappingEnhancers/nCommonVariantsTotal>0.001)
		max_enr = max(res_temp$enrichment)
		max_recall = max(res_temp$recall)
	}

	## plot enrichment heatmap
	plot_heatmap(res, "enrichment", "Enrichment", enr_ovl_colors, max_enr, out_enr, p_threshold)

	## plot recall heatmap
	plot_heatmap(res, "recall", "Recall", recall_ovl_colors, max_recall, out_recall, p_threshold)
}

main()
