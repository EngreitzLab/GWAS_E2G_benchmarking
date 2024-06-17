suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(ggplot2)
  library(ggpubr)})
		
result_files = snakemake@input$overlap_results %>% as.character() %>% strsplit(" ") %>% unlist()
trait_list = snakemake@params$trait_list %>% as.character() %>% strsplit(" ") %>% unlist()
trait_group = snakemake@wildcards$trait_group
cp = fread(snakemake@input$color_palette)
out_plot = snakemake@output$out_plot
out_table = snakemake@output$out_table

# read in data
for (i in 1:length(result_files)){
	temp = fread(result_files[i])
	temp = dplyr::filter(temp, trait %in% trait_list)
	if (i==1) {df = temp} else {
		df = rbind(df, temp)
	}
}

# aggregate, (mean for metrics, average variances 
df = mutate(df, var_enrichment = sd_enrichment^2, var_precision = sd_precision^2, var_recall = sd_recall^2) %>%
	group_by(method) %>%
	summarize(n_traits = n(), n_variants = sum(num_variants), enrichment = mean(enrichment), precision = mean(precision), recall=mean(recall),
						var_enrichment=mean(var_enrichment), var_precision=mean(var_precision), var_recall=mean(var_recall)) %>%
	mutate(sd_enrichment = sqrt(var_enrichment), sd_precision = sqrt(var_precision), sd_recall = sqrt(var_recall)) %>%
	left_join(cp, by="method")

# plot
pred_colors = cp$hex
names(pred_colors) = cp$pred_name_long
n_variants = as.numeric(df$n_variants) %>% mean()

enr_rec = ggplot(df, aes(x=recall, y=enrichment, xmin=recall-sd_recall, xmax=recall+sd_recall, ymin=enrichment-sd_enrichment, ymax=enrichment+sd_enrichment, color=pred_name_long)) +
	geom_pointrange() +
	scale_color_manual(values=pred_colors, name="Predictor") +
	labs(x="Recall", y="Enrichment") +
	#coord_fixed() +
	theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.text=element_text(size=7),
		legend.title=element_text(size=8), legend.position='None')

prec_rec = ggplot(df, aes(x=recall, y=precision, xmin=recall-sd_recall, xmax=recall+sd_recall, ymin=precision-sd_precision, ymax=precision+sd_precision, color=pred_name_long)) +
	geom_pointrange() +
	scale_color_manual(values=pred_colors, name="Predictor") +
	labs(x="Recall", y="Precision") +
	ylim(c(0, 1)) +
	#coord_fixed() +
	theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.text=element_text(size=7),
		legend.title=element_text(size=8), legend.position='None')

# aggregate plots and save things
g = ggarrange(enr_rec, prec_rec, nrow=1, common.legend=TRUE, legend="right")
g = annotate_figure(g, top=text_grob(paste0("Trait group: ", trait_group, " (N=", n_variants, ")"),size=10))
ggsave(out_plot, g, height=5, width=10)

fwrite(df, out_table, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
