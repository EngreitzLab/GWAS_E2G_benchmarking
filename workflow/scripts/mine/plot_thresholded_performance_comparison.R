suppressPackageStartupMessages({
		library(ggplot2)
        library(dplyr)
        library(tidyr)
    	library(tibble)
        library(cowplot)
		library(data.table)
	source("/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/workflow/scripts/helpful_math.R")})

filter_linking_data <- function(df, map, trait_groups, biosample_groups, type){
	filt_list = vector("list", nrow(map))
	if (type=="enhancers"){
		for (i in 1:nrow(map)){
			biosample_name = map$biosample[i]; trait_name = map$trait[i]
			message(biosample_name)
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
		legend.text=element_text(size=7), legend.title=element_text(size=8), legend.position='None', aspect.ratio = 1)

	return(enr_rec)

}

plot_gene_linking_one_set<- function(df_enhancers, df_baseline, withPoPS, p_threshold, pred_colors, out_table){
	# process predictor df
	df_input = df_enhancers
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

	prec_rec = ggplot(df, aes(x=recall, y=precision, color=pred_name_long)) +
		geom_point(data = df_other, alpha = 0.25, shape = 16, size = 3) +
		geom_pointrange(aes(xmin=CI_recall_low, xmax=CI_recall_high)) +
		geom_linerange(aes(ymin=CI_precision_low, ymax=CI_precision_high)) +
		scale_color_manual(values=pred_colors, name="Predictor") +
		labs(x=x_label, y="Precision (fraction of correct\npredicted credible set-gene links)") +
		xlim(c(0, x_max)) + ylim(c(0, y_max)) +
		ggtitle(label=plot_title) +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size=7), plot.title=element_text(size=10), legend.position = "None", aspect.ratio = 1)
	if (withPoPS) {
		prec_rec = prec_rec + theme(legend.text=element_text(size=7), legend.title=element_text(size=8), legend.position='right', aspect.ratio = 1)
	}
		

	return(prec_rec)
}


### MAIN
results_dir_all <- "/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results/2024_0908_sc_benchmarking"
methods_all <- c("scE2G_multiome", "scE2G_ATAC", "scABC", "ABC_distanceToTSS", "STARE", "SnapATAC", "Signac", "Cicero", "FigR", "SCENT", "ScenicPlus",
	 "DIRECTNET", "ArchR")
varOvl_files <- file.path(results_dir_all, methods_all, "variant_overlap",  "enrichmentRecall.thresholded.traitByBiosample.tsv.gz")
geneLinking_files <- file.path(results_dir_all, methods_all, "gene_linking", "precisionRecall.traitByBiosample.tsv.gz")

results_dir_hd <- "/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results/2024_0911_hd_blood"
varOvl_hd <- file.path(results_dir_hd, "scE2G_multiome", "variant_overlap",  "enrichmentRecall.thresholded.traitByBiosample.tsv.gz")
geneLinking_hd <-  file.path(results_dir_hd, "scE2G_multiome", "gene_linking", "precisionRecall.traitByBiosample.tsv.gz")

out_dir <- "/oak/stanford/groups/engreitz/Users/sheth/scE2G_analysis/2024_0823_GWAS_variant_interpretation/benchmarking"
out_plots <- file.path(out_dir, "thresholded_performance_comparison.scatter.pdf")
out_overlap <- file.path(out_dir, "thresholded_performance_comparison.variant_overlap.tsv")
out_linking <- file.path(out_dir, "thresholded_performance_comparison.gene_linking.tsv")

df_baseline <- fread(file.path(results_dir_all, "baseline", "gene_linking", "precisionRecall.byTrait.tsv.gz"))

cp <- fread(file.path(results_dir_all, "plots", "colorPalette.tsv"))
cp_add <- data.frame(method = "scE2G_multiome_hd", pred_name_long = "scE2G (Multiome, granular clusters)", hex = "#430b4e")
cp <- rbind(cp, cp_add)

p_threshold <- 0.05

# read in data for sc_benchmarking
trait_groups <- fread(file.path(results_dir_hd, "reference_configs", "trait_groups.tsv"))
biosample_groups <- fread(file.path(results_dir_hd, "reference_configs", "biosample_groups.tsv"))
biosample_names <- c("PBMC5", "BMMC5", "BMMC5_Erythroid")
trait_names <- c("Blood_counts", "Blood_counts", "RBC_traits")	

map = data.frame(biosample = biosample_names, trait = trait_names)
map = mutate(map, key = paste0(biosample, ".", trait))

df_overlap = lapply(varOvl_files, fread, sep="\t") %>% rbindlist() %>% as_tibble() %>%
	mutate(key = paste0(biosample, ".", trait)) %>%
	dplyr::filter(key %in% map$key) %>%
	left_join(cp, by="method")

df_linking = lapply(geneLinking_files, fread, sep="\t") %>% rbindlist() %>% as_tibble() %>%
	filter_linking_data(map, trait_groups, biosample_groups, "enhancers") %>%
	left_join(cp, by="method")

df_baseline = as_tibble(df_baseline) %>%
	filter_linking_data(map, trait_groups, biosample_groups, "baseline") %>%
	left_join(cp, by="method")

# add hd_blood data
map_hd <- fread("/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/config/config_comparisons_hd_blood.tsv") %>%
	dplyr::filter(name %in% c("Blood_counts_matched", "RBC_comparison")) %>%
	dplyr::select(-name) %>%
	mutate(key = paste0(biosample, ".", trait))

trait_groups_hd <- fread(file.path(results_dir_hd, "reference_configs", "trait_groups.tsv"))
biosample_groups_hd <- fread(file.path(results_dir_hd, "reference_configs", "biosample_groups.tsv"))

df_overlap_hd = fread(varOvl_hd)  %>% as_tibble() %>%
	mutate(key = paste0(biosample, ".", trait)) %>%
	dplyr::filter(key %in% map_hd$key) %>%
	mutate(method = "scE2G_multiome_hd") %>%
	left_join(cp, by="method")

df_linking_hd = fread(geneLinking_hd) %>% as_tibble() %>%
	filter_linking_data(map_hd, trait_groups_hd, biosample_groups_hd, "enhancers") %>%
	mutate(method = "scE2G_multiome_hd") %>%
	left_join(cp, by="method")

# combine
df_overlap <- rbind(df_overlap, df_overlap_hd)
df_linking <- rbind(df_linking, df_linking_hd)

## make panels
pred_colors = cp$hex
names(pred_colors) = cp$pred_name_long

overlap = plot_variant_overlap(df_overlap, p_threshold, pred_colors, out_overlap)
linking = plot_gene_linking_one_set(df_linking, df_baseline, FALSE, p_threshold, pred_colors, out_linking)
with_PoPS = plot_gene_linking_one_set(df_linking, df_baseline, TRUE, p_threshold, pred_colors, out_linking)
#linking = plot_gene_linking(df_linking, df_baseline, p_threshold, pred_colors, out_linking)

## assemble panels
grid = cowplot::plot_grid(overlap, linking, with_PoPS, nrow=1, ncol=3, align = "h", rel_widths=c(1,1,2.3))
ggsave(out_plots, grid, width=12, height=3.5)

