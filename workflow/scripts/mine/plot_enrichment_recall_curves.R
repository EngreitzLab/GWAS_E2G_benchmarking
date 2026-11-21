suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(colorspace)
	library(data.table)
	library(tidyr)
	library(cowplot)
	source("/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/workflow/scripts/helpful_math.R")})

add_plotting_params <- function(df, cp, p_threshold) {
	df = dplyr::filter(df, recall>0.005, p_adjust_enr<p_threshold, !is.na(enrichment),  # remove 0,0 points ( and arbitrarily low recall), non-significant and NA enrichments
		nVariantsTotal>25, nCommonVariantsOverlappingEnhancers/nCommonVariantsTotal>0.001) %>% # ...and cases with very few enhancers or variants
		left_join(cp, by=c("method")) # add hex and pred_name_long and score_threshold
		
	# add n_points
	n = dplyr::select(df, method, biosample, trait, threshold) %>%
		group_by(method, biosample, trait) %>%
		summarize(n_points = n())
	print(n)

	df = left_join(df, n, by=c("method", "biosample", "trait"))

	return(df)
}

plot_single_panel <- function(df, key_this, pred_colors) {
	# organize
	df = dplyr::filter(df, key==key_this)
	df_thresh = dplyr::filter(df, threshold==score_threshold) # suggested thresholds
	df_binary = dplyr::filter(df, n_points<=2, threshold!=0)
	df_plot = dplyr::filter(df, n_points>2)

	# figure out "N=""
	print(key_this)
	print(summary(df))
	n_label = get_n_equals(df, "variant_overlap")
	x_label = paste0("Recall (fraction of variants in predicted enhancers)\n", n_label, " fine-mapped variants")

	# axis limits
	max_recall =  ifelse(max(df$recall, na.rm=TRUE)<=0, 1, max(df$recall, na.rm=TRUE))
	abs_max_enr = max(df$enrichment, na.rm=TRUE)
	if (abs_max_enr<=0){
		max_enr = 10
	} else {
		max_enr_90 = unname(quantile(df$enrichment, probs=0.9, na.rm = TRUE))[1]
		if (abs_max_enr-max_enr_90>10){ # true max is an outlier
			max_enr = max_enr_90
		} else {
			max_enr = abs_max_enr
		}
	}

	g=ggplot(data=df_plot, aes(x=recall, y=enrichment, color=pred_name_long)) +
		geom_line(linewidth=1) +
		geom_point(data=df_binary, aes(x=recall, y=enrichment, color=pred_name_long), size=4) +
		geom_point(data=df_thresh, aes(x=recall, y=enrichment, color=pred_name_long), size=4) +
		geom_linerange(data=df_thresh, aes(ymin=CI_enr_low, ymax=CI_enr_high)) +
		scale_color_manual(values=pred_colors) +
		ylab("Enrichment (GWAS variants vs. common variants)") + xlab(x_label) +
		labs(col="Predictor") +
		xlim(c(0, max_recall)) + ylim(c(0, max_enr)) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), plot.title=element_text(size=10),
			legend.text = element_text(size=7), legend.title=element_text(size=8), legend.position="right", legend.direction="vertical",
			aspect.ratio = 1)
	
	return(g)
}

## input files + params
results_dir_all <- "/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results/2024_0908_sc_benchmarking"
results_dir_hd <- "/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results/2024_0911_hd_blood"
out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0823_GWAS_variant_interpretation/benchmarking"
out_combined <- file.path(out_dir, "combined_all_hd_enrichmentRecall.pdf")
out_singles <- file.path(out_dir, "singles_all_hd_enrichmentRecall.pdf")

ER_out_all <- c(file.path(results_dir_all, "plots", "enrichmentRecallCurves", "Blood_counts_coarse.values.tsv.gz"),
	file.path(results_dir_all, "plots", "enrichmentRecallCurves", "Erythroid_comparison.values.tsv.gz"))
ER_out_hd <- c(file.path(results_dir_hd, "plots", "enrichmentRecallCurves", "Blood_counts_matched.values.tsv.gz"),
	file.path(results_dir_hd, "plots", "enrichmentRecallCurves", "RBC_comparison.values.tsv.gz"))
p_threshold = 0.05

cp <- fread(file.path(results_dir_all, "plots", "colorPalette.tsv"))
cp_add <- data.frame(method = "scE2G_multiome_hd", pred_name_long = "scE2G (Multiome, granular clusters)", hex = "#430b4e")
cp <- rbind(cp, cp_add)

df_sep_all <- lapply(ER_out_all, fread) %>% rbindlist %>%
	dplyr::filter(key != "all_matched.all_matched")
df_sep_add <- lapply(ER_out_hd, fread) %>% rbindlist %>%
	dplyr::filter(key != "all_matched.all_matched") %>%
	mutate(method = "scE2G_multiome_hd")
df_sep <- rbind(df_sep_all, df_sep_add)

scores <- dplyr::select(df_sep, method, score_threshold) %>% distinct()
map <- dplyr::select(df_sep, biosample, trait, key) %>% distinct()
df_sep <- dplyr::select(df_sep, -c(score_threshold, hex, n_points, pred_name_long))

# ER_files = snakemake@input$ER %>% strsplit(" ")
# cp = fread(snakemake@input$colors, sep="\t") # method, pred_name_long, hex
# method_names = snakemake@params$method_names %>% strsplit(" ") %>% unlist()
# score_thresholds = snakemake@params$score_thresholds %>% as.character() %>% strsplit(" ") %>% unlist() %>% as.numeric()
# biosample_names = snakemake@params$biosample_names %>% strsplit(" ") %>% unlist()
# trait_names = snakemake@params$trait_names %>% strsplit(" ") %>% unlist()
# p_threshold = snakemake@params$p_threshold %>% as.numeric()
# comparison = snakemake@wildcards$comparison_name
# out_combined = snakemake@output$out_combined
# out_singles = snakemake@output$out_singles
# out_table = snakemake@output$out_table

# # read in tables
# map = data.frame(biosample=biosample_names, trait=trait_names)
# map = mutate(map, key = paste0(biosample, ".", trait))

# df_sep = lapply(ER_files, fread, sep="\t") %>%
# 	rbindlist() %>% as_tibble() %>%
# 	mutate(key = paste0(biosample, ".", trait)) %>%
# 	dplyr::filter(key %in% map$key) %>%
# 	dplyr::select(-c(biosampleGroup, traitGroup))

# aggregated metrics for this comparison
df_agg = group_by(df_sep, threshold, method)
df_agg = summarize_grouped_enrichment_recall(df_agg, p_threshold)
df_agg = mutate(df_agg, 
	trait = "all_matched", biosample = "all_matched", key=paste0(trait, ".", biosample))

# make plotting params df
#scores = data.frame(method=method_names, score_threshold=score_thresholds)
cp = left_join(cp, scores, by="method")

# format and save dfs
df_sep = add_plotting_params(df_sep, cp, p_threshold)
df_agg = add_plotting_params(df_agg, cp, p_threshold)
df_all = rbind(df_sep, df_agg)
#fwrite(df_all, out_table, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

## plotting
pred_colors = cp$hex
names(pred_colors) = cp$pred_name_long

# aggregate plot
a = plot_single_panel(df_agg, "all_matched.all_matched", pred_colors)
legend_v = as_grob(cowplot::get_plot_component(a, "guide-box-right"))
ggsave(out_combined, a, width=6, height=4)

# individual plots
plot_list = vector("list", nrow(map))
for (i in 1:nrow(map)){
	title = paste0("Biosample: ", map$biosample[i], "\nTrait: ", map$trait[i])
	plot_list[[i]] = plot_single_panel(df_sep, map$key[i], pred_colors) + ggtitle(title) + theme(legend.position="None")
}

# make grid...
if (nrow(map)==1){
	g = plot_list[[1]] + theme(legend.position="right")
	ggsave(out_singles, g, width=6, height=4)

} else {
	if (nrow(map)==2) {
		grid = plot_grid(plot_list[[1]], plot_list[[2]], legend_v, nrow=1, ncol=3)
		w = 9
		h = 4
	} else {
		top_row = plot_grid(plot_list[[1]], legend_v, nrow=1, ncol=2)
		n_plots = nrow(map)
		n_extra = nrow(map) - 1 # beyond row 1
		n_extra_rows = ceiling(n_extra/2)
		if (n_extra %% 2 > 0 ){
			plot_list[[nrow(map)+1]] = NULL
			n_plots = n_plots+1
		}
		etc = plot_grid(plotlist=plot_list[-1], nrow=n_extra_rows, ncol=2)
		grid = plot_grid(top_row, etc, nrow=2, ncol=1, rel_heights=c(1, n_extra_rows))
		w = 6
		h = 3*(n_extra_rows+1)
	}
	ggsave(out_singles, grid, width=w, height=h)
}
