import os
import polars as pl

# color palette 	
rule generate_color_palette:
	params:
		names = config["methods"],
		methods_config = config["methodsTable"]
	output:
		colorPalette = os.path.join(RESULTS_DIR, "plots", "colorPalette.tsv"),		
	conda: 
		os.path.join(ENV_DIR, "GWAS_env.yml")
	script: 
		os.path.join(SCRIPTS_DIR, "visualization", "color_palette.R")

# heatmaps
rule plot_enrichment_recall_heatmaps:
	input:
		variantOverlap = os.path.join(RESULTS_DIR, "{method}", "variant_overlap", "enrichmentRecall.traitByBiosample.tsv.gz"),
	params:
		p_threshold = config["thresholdPval"]
	output:
		outFile_enrichment =  os.path.join(RESULTS_DIR, "plots", "enrichmentHeatmaps", "{method}.enrichmentHeatmap.withMetrics.pdf"),
		outFile_recall = os.path.join(RESULTS_DIR, "plots", "recallHeatmaps", "{method}.recallHeatmap.withMetrics.pdf")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(ENV_DIR, "GWAS_env.yml")
	script:
		os.path.join(SCRIPTS_DIR, "visualization", "plot_enrichment_recall_heatmaps.R")



# rule plot_performance: # per comparison
# 	input:
# 		overlaps = [os.path.join(config["outDir"], method, "tissues", "{tissue}", "{trait_group}_overlap_results.tsv") for method in config["methods"]]
# 	params:
# 		methods = config["methods"],
# 		methodsConfig = config["methodsTable"] # to get colors
# 	conda: 
# 		os.path.join(config["envDir"], "GWAS_env.yml")
# 	resources:
# 		mem_mb = determine_mem_mb
# 	output:
# 		outPlot = os.path.join(config["outDir"], "plots",  "Tissue{tissue}_Trait{trait_group}_overlap_results.pdf"),
# 		outFile = os.path.join(config["outDir"], "plots",  "Tissue{tissue}_Trait{trait_group}_overlap_results.tsv")
# 	script:
# 		os.path.join(config["codeDir"], "plot_performance_metrics.R")