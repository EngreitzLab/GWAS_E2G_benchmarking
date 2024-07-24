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

		print(filt_list[[i]])
	}

	filt = rbindlist(filt_list) %>% as_tibble() %>% distinct()
	return(filt)
}

plot_variant_overlap <- function(df, p_threshold, pred_colors, out_table){
	df = group_by(df, method, pred_name_long, hex) %>%
		summarize_grouped_enrichment_recall(p_threshold)
	fwrite(df, out_table, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)

	# get N = for label
	n_label = get_n_equals(df, "variant_overlap")
	x_label = paste0("Recall (fraction of variants in predicted enhancers) \n", n_label, " fine-mapped variants")

	# plot!
	df = mutate(df,
		CI_enr_low = ifelse(CI_enr_low<0, 0, CI_enr_low),
		CI_recall_low = ifelse(CI_recall_low<0, 0, CI_recall_low))

	enr_rec = ggplot(df, aes(x=recall, y=enrichment, color=pred_name_long)) +
		geom_pointrange(aes(xmin=CI_recall_low, xmax=CI_recall_high)) +
		geom_pointrange(aes(ymin=CI_enr_low, ymax=CI_enr_high)) +
		scale_color_manual(values=pred_colors, name="Predictor") +
		labs(x=x_label, y="Enrichment\n(GWAS variants vs. common variants)") +
		xlim(c(0, max(df$CI_recall_high))) + ylim(c(0, max(df$CI_enr_high))) +
		ggtitle(label="Variant overlap") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 7), plot.title=element_text(size=10),
		legend.text=element_text(size=7), legend.title=element_text(size=8), legend.position='None')

	return(enr_rec)

}

# get number of credible sets per trait and a data table of unique credible sets
plot_gene_linking <- function(df_enhancers, df_baseline, p_threshold, pred_colors, out_table){
	# process predictor df
	df_enhancers = group_by(df_enhancers, method, pred_name_long, intersectPoPS) %>%
		summarize_grouped_precision_recall(p_threshold, "enhancers") %>%
		dplyr::select(-c(nCredibleSetsOverlappingEnhancersCorrectGene, nCredibleSetsOverlappingEnhancersAnyGene))

	# process baseline df
	df_baseline = group_by(df_baseline, method, pred_name_long) %>%
		summarize_grouped_precision_recall(p_threshold, "baseline") %>%
		dplyr::select(-c(nCredibleSetGenePredictionsCorrect, nCredibleSetGenePredictions))

	# combine + save
	df = rbind(df_enhancers, df_baseline)
	fwrite(df, out_table, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)

	# make paired df for lines..
	df_pairs = dplyr::select(df_enhancers, method, pred_name_long, precision, recall, intersectPoPS) %>%
		mutate(intersectPoPS = ifelse(intersectPoPS, "withPoPS", "withoutPoPS")) %>%
		pivot_wider(names_from=intersectPoPS, values_from=c(precision, recall), names_sep="_") %>%
		dplyr::filter(precision_withPoPS!=precision_withoutPoPS || recall_withPoPS != recall_withPoPS)
	print(df_pairs)

	# get N = for label
	n_label = get_n_equals(df, "gene_linking")
	x_label = paste0("Recall (fraction of credible sets linked to target gene)\n", n_label, " credible sets")

	# plot!
	df = mutate(df,
		CI_precision_low = ifelse(CI_precision_low<0, 0, CI_precision_low),
		CI_recall_low = ifelse(CI_recall_low<0, 0, CI_recall_low),
		CI_precision_high = ifelse(CI_precision_high>1, 1, CI_precision_high))

	prec_rec = ggplot(df, aes(x=recall, y=precision, color=pred_name_long)) +
		geom_curve(data=df_pairs, aes(x=recall_withoutPoPS, xend=recall_withPoPS, y=precision_withoutPoPS, yend=precision_withPoPS, color=pred_name_long),
			linetype="82", linewidth=0.35, alpha=1, curvature=0.55, lineend="butt", show.legend=FALSE,
			arrow=arrow(length=unit(0.075, "npc"), angle=30, type="closed")) +
		geom_pointrange(aes(xmin=CI_recall_low, xmax=CI_recall_high)) +
		geom_pointrange(aes(ymin=CI_precision_low, ymax=CI_precision_high)) +
		scale_color_manual(values=pred_colors, name="Predictor") +
		labs(x=x_label, y="Precision (fraction of correct\npredicted credible set-gene links)") +
		xlim(c(0, max(df$CI_recall_high))) + ylim(c(0, max(df$CI_precision_high))) +
		ggtitle(label="Linking variants to genes") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size=7), plot.title=element_text(size=10),
		legend.text=element_text(size=7), legend.title=element_text(size=8), legend.position='right')

	return(prec_rec)


}

### MAIN
varOvl_files = snakemake@input$variantOverlap %>% strsplit(" ")
geneLinking_files = snakemake@input$geneLinking %>% strsplit(" ")
df_baseline = fread(snakemake@input$geneLinking_baseline, sep="\t")
cp = fread(snakemake@input$colors, sep="\t")
trait_groups = fread(snakemake@input$traitGroups, sep="\t")
biosample_groups = fread(snakemake@input$biosampleGroups, sep="\t")
biosample_names = snakemake@params$biosample_names %>% strsplit(" " ) %>% unlist()
trait_names = snakemake@params$trait_names %>% strsplit(" ") %>% unlist()
p_threshold = snakemake@params$p_threshold %>% as.numeric()
out_plots = snakemake@output$plots
out_overlap =snakemake@output$table_overlap
out_linking = snakemake@output$table_linking

## read in data
map = data.frame(biosample = biosample_names, trait = trait_names)
map = mutate(map, key = paste0(biosample, ".", trait))

df_overlap = lapply(varOvl_files, fread, sep="\t") %>% rbindlist() %>% as_tibble() %>%
	mutate(key = paste0(biosample, ".", trait)) %>%
	dplyr::filter(key %in% map$key) %>%
	left_join(cp, by="method")
print(df_overlap)

df_linking = lapply(geneLinking_files, fread, sep="\t") %>% rbindlist() %>% as_tibble() %>%
	filter_linking_data(map, trait_groups, biosample_groups, "enhancers") %>%
	left_join(cp, by="method")
print(df_linking)

df_baseline = as_tibble(df_baseline) %>%
	filter_linking_data(map, trait_groups, biosample_groups, "baseline") %>%
	left_join(cp, by="method")

## make panels
pred_colors = cp$hex
names(pred_colors) = cp$pred_name_long

overlap = plot_variant_overlap(df_overlap, p_threshold, pred_colors, out_overlap)
linking = plot_gene_linking(df_linking, df_baseline, p_threshold, pred_colors, out_linking)

## assemble panels
grid = plot_grid(overlap, linking, nrow=1, ncol=2, rel_widths=c(1,2))
ggsave(out_plots, grid, width=8, height=3)

