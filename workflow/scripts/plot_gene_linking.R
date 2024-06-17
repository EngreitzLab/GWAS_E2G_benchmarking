suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(ggplot2)
  library(ggpubr)})
		
result_files = snakemake@input$overlap_results %>% as.character() %>% strsplit(" ") %>% unlist()
trait_group = snakemake@wildcards$trait_group
cp = fread(snakemake@input$color_palette)
out_plot = snakemake@output$out_plot
out_table = snakemake@output$out_table

# read in data
for (i in 1:length(result_files)){
	temp = fread(result_files[i])
	if (i==1) {df = temp} else {
		df = rbind(df, temp)
	}
}

# format for plotting (colors!)
df = left_join(df, cp, by="method")
df$num_pops_genes[df$use_pops!="withPoPS"] = 0
df$key_vals = paste0("(# PoPS genes = ", df$num_pops_genes, ", # pred. genes = ", df$num_pred_genes, ")")
df$key = paste0(df$pred_name_long, "\n", df$key_vals)

cp_new = dplyr::select(df, method, key, hex) %>% distinct() 
methods = unique(cp_new$method)
if (length(methods) < length(unique(cp_new$key))) {
	# set new hex codes
	for (i in 1:length(methods)){
		cp_this =  dplyr::filter(cp_new, method==methods[i])
		n = nrow(cp_this)
		if (n>1){
			hex_base = cp_this$hex[1]
			cols = colorRampPalette(c(hex_base, '#ffffff'))(n+2) # range of colors from hex_base to white (where white is 2 above the base)
			for (k in 1:nrow(cp_this)){
				cp_new$hex[cp_new$key==cp_this$key[k]] = cols[k]  # set corresponding hex code
			}
		}
	}
}



# plot
pred_colors = cp_new$hex
names(pred_colors) = cp_new$key
n_variants = as.numeric(df$num_relevant) %>% mean()

# enr_rec = ggplot(df, aes(x=recall, y=enrichment, xmin=recall-sd_recall, xmax=recall+sd_recall, ymin=enrichment-sd_enrichment, ymax=enrichment+sd_enrichment, color=key)) +
# 	geom_pointrange() +
# 	scale_color_manual(values=pred_colors, name="Predictor") +
# 	labs(x="Recall", y="Enrichment") +
# 	#coord_fixed() +
# 	theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.text=element_text(size=7),
# 		legend.title=element_text(size=8), legend.position='None')

prec_rec = ggplot(df, aes(x=recall, y=precision, xmin=recall-sd_recall, xmax=recall+sd_recall, ymin=precision-sd_precision, ymax=precision+sd_precision, color=key)) +
	geom_pointrange() +
	scale_color_manual(values=pred_colors, name="Predictor") +
	labs(x="Recall", y="Precision") +
	ylim(c(0, 1)) +
	#coord_fixed() +
	theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.text=element_text(size=7),
		legend.title=element_text(size=8), legend.position='None')

# aggregate plots and save things
#g = ggarrange(enr_rec, prec_rec, nrow=1, common.legend=TRUE, legend="bottom")
g = ggarrange(NULL, prec_rec, nrow=1, common.legend=TRUE, legend="right", legend.direction="vertical")
g = annotate_figure(g, top=text_grob(paste0("Trait group: ", trait_group, " (N=", n_variants, ")"),size=10))
ggsave(out_plot, g, height=5, width=12)

df = dplyr::select(df, -key) # remove this column bc new line
fwrite(df, out_table, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
