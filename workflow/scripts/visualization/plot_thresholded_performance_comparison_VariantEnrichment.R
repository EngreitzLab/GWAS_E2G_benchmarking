suppressPackageStartupMessages({
		library(ggplot2)
        library(dplyr)
        library(tidyr)
    	library(tibble)
        library(cowplot)
		library(data.table)
		source(snakemake@params$helpful_math)})

filter_linking_data <- function(df, map, trait_groups, biosample_groups, type){
	filt_list = vector("list", nrow(map))
	print(map)
	if (type=="enhancers"){
		for (i in 1:nrow(map)){
			biosample_name = map$biosample[i]; trait_name = map$trait[i]
		
			if (biosample_name=="ALL") {
				df_biosample = df
			} else if (biosample_name %in% colnames(biosample_groups)){
				df_biosample = dplyr::filter(df, biosample %in% biosample_groups[[biosample_name]])
			} else {
				df_biosample = dplyr::filter(df, biosample==biosample_name)
			}

			if (trait_name %in% colnames(trait_groups)){
				df_trait = dplyr::filter(df_biosample, trait %in% trait_groups[[trait_name]])
			} else {
				df_trait = dplyr::filter(df_biosample, trait==trait_name)
			}
			filt_list[[i]] = df_trait
		}
	} else if (type=="baseline"){
		for (i in 1:nrow(map)){
			trait_name = map$trait[i]
			if (trait_name %in% colnames(trait_groups)){
				df_trait = dplyr::filter(df, trait %in% trait_groups[[trait_name]])
			} else {
				df_trait = dplyr::filter(df, trait==trait_name)
			}
			filt_list[[i]] = df_trait
		}
	}

	filt = rbindlist(filt_list) %>% as_tibble() %>% distinct()
	return(filt)
}

plot_variant_overlap <- function(df, p_threshold, pred_colors, out_table, all_var, trait_groups){
	df_input <- df
	df = group_by(df, method, pred_name_long, hex) %>%
		summarize_grouped_enrichment_recall(p_threshold)
	fwrite(df, out_table, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)

	# get N = for label
	traits <- unique(df_input$trait)
	n_label = get_n_equals_variant_overlap(df_input, all_var, traits, trait_groups)
	x_label = paste0("Recall (fraction of variants in predicted enhancers) \n", n_label)

	# find max for y axis
	enr_sort <- sort(c(df$CI_enr_high, df$enrichment), decreasing = TRUE)
	if (enr_sort[1] > 1.2 * enr_sort[2]){ #is the max more than 20% of the second highest value?
		enr_max <- enr_sort[2] * 1.05
	} else {
		enr_max <- enr_sort[1]
	}

	# plot!
	df = mutate(df,
		CI_enr_low = ifelse(CI_enr_low<0, 0, CI_enr_low),
		CI_recall_low = ifelse(CI_recall_low<0, 0, CI_recall_low))

	enr_rec = ggplot(df, aes(x=recall, y=enrichment, color=pred_name_long)) +
		geom_linerange(aes(xmin=CI_recall_low, xmax=CI_recall_high)) +
		geom_linerange(aes(ymin=CI_enr_low, ymax=CI_enr_high)) +
		geom_point(shape = 16, size = 4) +
		scale_color_manual(values=pred_colors, name="Predictor") +
		labs(x=x_label, y="Enrichment\n(GWAS variants vs. common variants)") +
		xlim(c(0, max(df$CI_recall_high))) + ylim(c(0, enr_max)) +
		ggtitle(label="Variant overlap") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 7), plot.title=element_text(size=10),
		legend.text=element_text(size=7), legend.title=element_text(size=8), legend.position='None', aspect.ratio = 1)

	return(enr_rec)

}


plot_gene_linking_one_set<- function(df_enhancers, df_baseline, withPoPS, p_threshold, pred_colors, out_table){
	df_input <- df_enhancers

	# process predictor df
	df_enhancers = group_by(df_enhancers, method, pred_name_long, intersectPoPS) %>%
		summarize_grouped_precision_recall(p_threshold, "enhancers") %>%
		dplyr::select(-c(nCredibleSetsOverlappingEnhancersCorrectGene, nCredibleSetsOverlappingEnhancersAnyGene))

	# process baseline df
	df_baseline = group_by(df_baseline, method, pred_name_long) %>%
		summarize_grouped_precision_recall(p_threshold, "baseline") %>%
		dplyr::select(-c(nCredibleSetGenePredictionsCorrect, nCredibleSetGenePredictions))

	# combine + save
	if (withPoPS){
		df_save = rbind(df_enhancers, df_baseline)
		fwrite(df_save, out_table, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
	}

	# get axis limits (constant for plus and minus PoPS)
	x_max <- min(1, max(df_enhancers$CI_recall_high))
	y_max <- min(1, max(df_enhancers$CI_precision_high))

	# pull relevant data for this plot
	df_enhancers = mutate(df_enhancers,
		CI_precision_low = ifelse(CI_precision_low<0, 0, CI_precision_low),
		CI_recall_low = ifelse(CI_recall_low<0, 0, CI_recall_low),
		CI_precision_high = ifelse(CI_precision_high>1, 1, CI_precision_high))

	df = dplyr::filter(df_enhancers, intersectPoPS == withPoPS)
	df_other <- dplyr::filter(df_enhancers, intersectPoPS == !withPoPS)

	# get N = for label #n_label = get_n_equals(df, "gene_linking")
	n_label <- get_n_equals_gene_linking(df_input)
	x_label = paste0("Recall (fraction of credible sets linked to target gene)\n", n_label)
	plot_title = ifelse(withPoPS, "Linking variants to known genes\n(E2G method int. PoPS)", "Linking variants to known genes\n(E2G method)")

	# plot!

	prec_rec = ggplot(df, aes(x=recall, y=precision, color=pred_name_long)) +
		geom_point(data = df_other, alpha = 0.25, shape = 16, size = 4) +
		geom_point(shape = 16, size = 4) + 
		geom_linerange(aes(xmin=CI_recall_low, xmax=CI_recall_high)) +
		geom_linerange(aes(ymin=CI_precision_low, ymax=CI_precision_high)) +
		scale_color_manual(values=pred_colors, name="Predictor") +
		labs(x=x_label, y="Precision (fraction of correct\npredicted credible set-gene links)") +
		xlim(c(0, x_max)) + ylim(c(0, y_max)) +
		ggtitle(label=plot_title) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size=7),
			plot.title=element_text(size=10), legend.position = "None", aspect.ratio = 1)
	if (withPoPS) {
		prec_rec = prec_rec + theme(legend.text=element_text(size=7), legend.title=element_text(size=8), legend.position='right', aspect.ratio = 1)
	}
		

	return(prec_rec)
}


### MAIN
varOvl_files = snakemake@input$variantOverlap %>% strsplit(" ")
# geneLinking_files = snakemake@input$geneLinking %>% strsplit(" ")
# df_baseline = fread(snakemake@input$geneLinking_baseline, sep="\t")
all_var = fread(snakemake@input$allVariants)
cp = fread(snakemake@input$colors, sep="\t")
trait_groups = fread(snakemake@input$traitGroups, sep="\t")
biosample_groups = fread(snakemake@input$biosampleGroups, sep="\t")
biosample_names = snakemake@params$biosample_names %>% strsplit(" " ) %>% unlist()
trait_names = snakemake@params$trait_names %>% strsplit(" ") %>% unlist()
p_threshold = snakemake@params$p_threshold %>% as.numeric()
out_plots = snakemake@output$plots
out_overlap =snakemake@output$table_overlap
# out_linking = snakemake@output$table_linking

## read in data
map = data.frame(biosample = biosample_names, trait = trait_names)
map = mutate(map, key = paste0(biosample, ".", trait))

df_overlap = lapply(varOvl_files, fread, sep="\t") %>% rbindlist() %>% as_tibble() %>%
	mutate(key = paste0(biosample, ".", trait)) %>%
	dplyr::filter(key %in% map$key) %>%
	left_join(cp, by="method")
print(df_overlap)

# df_linking = lapply(geneLinking_files, fread, sep="\t") %>% rbindlist() %>% as_tibble() %>%
# 	filter_linking_data(map, trait_groups, biosample_groups, "enhancers") %>%
# 	left_join(cp, by="method")
# print(df_linking)

# df_baseline = as_tibble(df_baseline) %>%
# 	filter_linking_data(map, trait_groups, biosample_groups, "baseline") %>%
# 	left_join(cp, by="method")

## make panels
pred_colors = cp$hex
names(pred_colors) = cp$pred_name_long

overlap = plot_variant_overlap(df_overlap, p_threshold, pred_colors, out_overlap, all_var, trait_groups)
# linking = plot_gene_linking_one_set(df_linking, df_baseline, FALSE, p_threshold, pred_colors, out_linking)
# with_PoPS = plot_gene_linking_one_set(df_linking, df_baseline, TRUE, p_threshold, pred_colors, out_linking)

# legend_v = as_grob(cowplot::get_plot_component(with_PoPS, "guide-box-right"))
# with_PoPS <- with_PoPS + theme(legend.position = "None")

## assemble panels
# grid = cowplot::plot_grid(overlap, linking, with_PoPS, nrow=1, ncol=3, align = "h", axis="tb", rel_widths=c(1,1,1))
# grid_l = plot_grid(grid, legend_v, nrow=1, ncol = 2, rel_widths=c(9,2))
# ggsave(out_plots, grid_l, width=10, height=3)

ggsave(out_plots, overlap, width = 10, height = 3)

