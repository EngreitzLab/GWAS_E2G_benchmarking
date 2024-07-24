suppressPackageStartupMessages({
		library(ggplot2)
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
	if(sum(sum(trait_dist))>0){
		order_traits = hclust(trait_dist, method = "ward.D2")$order
		traits_ordered = colnames(M)[order_traits]
	} else {
		traits_ordered = colnames(M)
	}


	biosample_dist = dist(1-cor(t(M)))
	biosample_dist[is.na(biosample_dist)] = 0
	if (sum(sum(biosample_dist))>0){
		order_biosamples = hclust(biosample_dist, method="ward.D2")$order
		biosamples_ordered = rownames(M)[order_biosamples]
	} else {
		biosamples_ordered = rownames(M)
	}


	res$trait = factor(res$trait, levels=traits_ordered, ordered=TRUE)
	res$biosample = factor(res$biosample, levels=biosamples_ordered, ordered=TRUE)
	return(res)
}

plot_heatmap <- function(res, metric, metric_name, method, colors, max_value, out_file){
	## subset data
	res$metric = res[[metric]]
	if (method=="intPoPS"){
		res = dplyr::filter(res, intersectPoPS)
	} else {
		res = dplyr::filter(res, !intersectPoPS)
	}

	## cluster
	res = cluster_traits_biosamples(res)

	## plotting params
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

	# cs_count
	n_cs = dplyr::select(res, trait, nCredibleSetsTotal) %>% distinct()
	cs_count = ggplot(n_cs, aes(x=trait, y=nCredibleSetsTotal)) +
		geom_bar(stat="identity", width=0.5, fill=bar_colors[1]) +
		ylab("# credible sets\nper trait") + xlab("") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.text.x = element_text(angle=60, hjust=1))

	# combine with cs count
	assembled = egg::ggarrange(hm, cs_count, nrow=2, ncol=1, heights=c(2, 0.2))
	ggsave(out_file, assembled, width=8, height=ht)	
}

main <- function() {
	## process data 
	res_file = snakemake@input$geneLinking
	ranges_file = snakemake@input$ranges
	p_threshold = snakemake@params$p_threshold %>% as.numeric()
	fixed_scale = snakemake@params$fixedScale
	out_precision = snakemake@output$outFile_precision
	out_recall = snakemake@output$outFile_recall
	out_precision_PoPS = snakemake@output$outFile_precision_PoPS
	out_recall_PoPS = snakemake@output$outFile_recall_PoPS

	enr_ovl_colors =  c("#f1eef6","#bdc9e1", "#74a9cf","#2b8cbe", "#045a8d") # PuBu
	recall_ovl_colors = c("#edf8fb","#b3cde3", "#8c96c6", "#8856a7", "#810f7c") # BuPu
	precision_colors = c("#edf8fb","#b2e2e2", "#66c2a4", "#2ca25f", "#006d2c")  # BuGn
	recall_link_colors = c("#ffffcc","#c2e699", "#78c679", "#31a354", "#006837") # YlGn
	

	res = fread(res_file, sep="\t", header=TRUE)
	res = dplyr::filter(res, is.finite(precision), is.finite(recall), !is.na(precision), !is.na(recall), nCredibleSetsTotal>0) # can filter to n>X to reduce noise maybe

	# get max values for color scale
	if (fixed_scale){
		ranges = fread(ranges_file, sep="\t")
		precision = dplyr::filter(ranges, metric=="precision_linking", min_max=="max"); max_precision = max(precision$value)
		precision_PoPS = dplyr::filter(ranges, metric=="precision_linking_intersectPoPS", min_max=="max"); max_precision_PoPS = max(precision_PoPS$value, na.rm=TRUE)
		recall = dplyr::filter(ranges, metric=="recall_linking", min_max=="max"); max_recall = max(recall$value) 
		recall_PoPS = dplyr::filter(ranges, metric=="recall_linking_intersectPoPS", min_max=="max"); max_recall_PoPS = max(recall_PoPS$value, na.rm=TRUE)
	} else {
		noPoPS = dplyr::filter(res, !intersectPoPS)
		max_precision = max(noPoPS$precision)
		max_recall = max(noPoPS$recall)

		withPoPS  = dplyr::filter(res, intersectPoPS)
		max_precision_PoPS = max(withPoPS$precision, na.rm=TRUE)
		max_recall_PoPS = max(withPoPS$recall, na.rm=TRUE)
	}
	


	# plot heatmaps
	for (metric in c("precision", "recall")){
		for (method in c("noPoPS", "intPoPS")){
				# choose output file
				if (metric=="precision"){
					colors = precision_colors
					metric_name = "Precision"
					if (method=="intPoPS"){
						out_file = out_precision_PoPS
						max_value = max_precision_PoPS
					} else {
						out_file = out_precision
						max_value = max_precision
					}
				} else {
					colors=recall_link_colors
					metric_name = "Recall"
					if (method=="intPoPS"){
						out_file = out_recall_PoPS
						max_value = max_recall_PoPS
					} else {
						out_file = out_recall
						max_value = max_recall
					}
				}

				plot_heatmap(res, metric, metric_name, method, colors, max_value, out_file)
		}
	}
}

main()
