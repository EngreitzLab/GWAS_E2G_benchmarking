suppressPackageStartupMessages({
		library(ggplot2)
        library(dplyr)
        library(scales)
        library(tidyr)
    	library(tibble)
        library(egg)
		library(data.table)})

### generate plot across predictors - for a given subset of biosample/traits 

## figure 2e: variant overlap enrichment / recall 
# K562 x MCV, HbA1c, RBC
# GM12878 x Lym, AID_Combined
plot_figure_2e <- function(df_ovl, p_threshold, pred_colors, out_plot, out_table){
	K562_traits =  c("MCV, HbA1c, RBC")
	GM12878_traits = c("Lym", "AID_Combined")

	# filter data 
	df_K562 = dplyr::filter(df_ovl, biosample== "K562", trait =="RBC_traits")
	#print(df_K562)

	df_GM12878 = dplyr::filter(df_ovl, biosample=="GM12878", trait=="WBC_traits")
	#print(df_GM12878)

	df = rbind(df_K562, df_GM12878)

	# get N for plotting label
	traits =  dplyr::select(df, trait, nVariantsTotal) %>% distinct()
	n_variants = sum(traits$nVariantsTotal)

	# aggregate data
	z = qnorm(p_threshold/2, lower.tail=FALSE) 
	df = mutate(df, SE_log_enr_sq = SE_log_enr^2, SE_recall_sq = SE_recall^2) %>%
		group_by(method) %>%
		summarize(enrichment = mean(enrichment),
			recall = mean(recall),
			SE_log_enr = sqrt(mean(SE_log_enr_sq)),
			SE_recall = sqrt(mean(SE_recall_sq)),
			nVariantsTotal = sum(nVariantsTotal)) %>%
		mutate(CI_enr_low = exp(log(enrichment)-z*SE_log_enr),
			CI_enr_high= exp(log(enrichment)+z*SE_log_enr),
			CI_recall_low = recall - z*SE_recall,
			CI_recall_high = recall + z*SE_recall)
	fwrite(df, out_table, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)

	# plot!
	df$CI_recall_low = pmax(df$CI_recall_low, 0)
	df$CI_enr_low = pmax(df$CI_enr_low, 0)
	enr_rec = ggplot(df, aes(x=recall, y=enrichment, color=method)) +
		geom_pointrange(aes(xmin=CI_recall_low, xmax=CI_recall_high)) +
		geom_pointrange(aes(ymin=CI_enr_low, ymax=CI_enr_high)) +
		scale_color_manual(values=pred_colors, name="Predictor") +
		labs(x=paste0("Recall (95% CI, N=", n_variants, " variants)"), y="Enrichment vs. common variants (95% CI)") +
		xlim(c(0, max(df$CI_recall_high))) + ylim(c(0, max(df$CI_enr_high))) +
		ggtitle(label="GWAS variant overlap", subtitle="K562 x MCV, HbA1c, RBC\nGM12878 x Lym, AID_Combined") +
		theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.text=element_text(size=7),
			legend.title=element_text(size=8), legend.position='right')

	ggsave(out_plot, enr_rec, width=4, height=4)
}

## figure 2f: gene-linking for 11 blood traits in blood biosamples
plot_figure_2f <- function(df_link, p_threshold, pred_colors, out_plot, out_table){

}

### MAIN
init_results_dir = "/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results/2024_0624_initial_test"
sc_results_dir = "/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results/2024_0630_sc_test"
methods = c("scE2G_multiome", "scE2G_ATAC", "SCENT", "ENCODE_rE2Gext", "ABC_full", "distanceToTSS", "PoPS")
pred_colors = c("#792374", "#006479", "#a5a083", "#9B241C", "#00488D", "#6e788d", "#1c2a43" )
names(pred_colors) = methods
p_threshold=0.05

# concatenate variant overlap results
ovl_methods = methods[1:5]
for (i in 1:length(ovl_methods)){
	file_name = file.path(sc_results_dir, methods[i], "variant_overlap",  "enrichmentRecall.traitByBiosample.tsv.gz")
	temp = fread(file_name, sep="\t")

	if (i==1){df_ovl = temp} else {df_ovl = rbind(df_ovl, temp)}
}

#n_traits = dplyr::select(df_ovl, trait, nVariantsTotal) %>% distinct()
#print(n_traits)

old_ovl = fread(file.path(init_results_dir, "ABC_full", "variant_overlap", "enrichmentRecall.traitByBiosample.tsv.gz"), sep="\t")
n_traits = dplyr::select(old_ovl, trait, nVariantsTotal) %>% distinct()
print(n_traits)

out_2e = file.path(sc_results_dir, "plots", "Fig2e_attempt_v2.pdf")
out_2e_tsv = file.path(sc_results_dir, "plots", "Fig2e_attempt_v2.tsv")
plot_figure_2e(df_ovl, p_threshold, pred_colors, out_2e, out_2e_tsv)