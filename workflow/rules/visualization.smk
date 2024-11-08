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

# gather range of each metric across predictors to define scales
rule gather_metric_ranges:
	input:
		variantOverlap = expand(os.path.join(RESULTS_DIR, "{method}", "variant_overlap", "enrichmentRecall.thresholded.traitByBiosample.tsv.gz"), method=config["methods"]),
		geneLinking =  expand(os.path.join(RESULTS_DIR, "{method}", "gene_linking", "precisionRecall.traitByBiosample.tsv.gz"), method=config["methods"])
	params:
		p_threshold = config["thresholdPval"]
	output:
		ranges = os.path.join(RESULTS_DIR, "plots", "metricRanges.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(ENV_DIR, "GWAS_env.yml")
	script:
		os.path.join(SCRIPTS_DIR, "visualization", "gather_metric_ranges.R")

# heatmaps for variant overlap enrichment + recall
rule plot_variant_overlap_heatmaps:
	input:
		variantOverlap = os.path.join(RESULTS_DIR, "{method}", "variant_overlap", "enrichmentRecall.thresholded.traitByBiosample.tsv.gz"),
		ranges = os.path.join(RESULTS_DIR, "plots", "metricRanges.tsv")
	params:
		p_threshold = config["thresholdPval"],
		fixedScale = config["plotFixedScale"],
	output:
		outFile_enrichment =  os.path.join(RESULTS_DIR, "plots", "enrichmentOverlapHeatmaps", "{method}.enrichmentOverlapHeatmap.withMetrics.pdf"),
		outFile_recall = os.path.join(RESULTS_DIR, "plots", "recallOverlapHeatmaps", "{method}.recallOverlapHeatmap.withMetrics.pdf")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(ENV_DIR, "GWAS_env.yml")
	script:
		os.path.join(SCRIPTS_DIR, "visualization", "plot_enrichment_recall_heatmaps.R")

# heatmaps for gene linking precision + recall
rule plot_gene_linking_heatmaps:
	input:
		geneLinking = os.path.join(RESULTS_DIR, "{method}", "gene_linking", "precisionRecall.traitByBiosample.tsv.gz"),
		ranges = os.path.join(RESULTS_DIR, "plots", "metricRanges.tsv")
	params:
		p_threshold = config["thresholdPval"],
		fixedScale = config["plotFixedScale"]
	output:
		outFile_precision =  os.path.join(RESULTS_DIR, "plots", "precisionLinkingHeatmaps", "{method}.precisionLinkingHeatmap.withMetrics.pdf"),
		outFile_recall = os.path.join(RESULTS_DIR, "plots", "recallLinkingHeatmaps", "{method}.recallLinkingHeatmap.withMetrics.pdf"),
		outFile_precision_PoPS =  os.path.join(RESULTS_DIR, "plots", "intPoPS.precisionLinkingHeatmaps", "{method}.intPoPS.precisionLinkingHeatmap.withMetrics.pdf"),
		outFile_recall_PoPS = os.path.join(RESULTS_DIR, "plots", "intPoPS.recallLinkingHeatmaps", "{method}.intPoPS.recallLinkingHeatmap.withMetrics.pdf"),
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(ENV_DIR, "GWAS_env.yml")
	script:
		os.path.join(SCRIPTS_DIR, "visualization", "plot_precision_recall_heatmaps.R")

# plot enrichment-recall curves (per comparison "name")
rule plot_enrichment_recall_curves:
	input: 
		# enr recall for each method for necessary biosamples
		ER = lambda wildcards: flatten([[os.path.join(RESULTS_DIR, method, "variant_overlap", f"enrichmentRecallAcrossThresholds.forBiosample{biosample}.tsv.gz")
			for biosample in get_biosample_names_per_method(comparisons.filter(pl.col("name")==wildcards.comparison_name), method)]
			for method in config["methods"]]),
		colors = os.path.join(RESULTS_DIR, "plots", "colorPalette.tsv"),
		allVariants = (os.path.join(RESULTS_DIR, "variants", "filteredGWASVariants.merged.sorted.tsv.gz")),
		traitGroups = os.path.join(RESULTS_DIR, "reference_configs", "trait_groups.tsv")
	params:
		method_names = config["methods"],
		biosample_names = lambda wildcards: comparisons.filter(pl.col("name")==wildcards.comparison_name)["biosample"].to_list(),
		trait_names = lambda wildcards: comparisons.filter(pl.col("name")==wildcards.comparison_name)["trait"].to_list(),
		score_thresholds = [get_col2_from_col1(methods_config, "method", method_name, "threshold") for method_name in config["methods"]],
		p_threshold = config["thresholdPval"],
		helpful_math = os.path.join(SCRIPTS_DIR, "helpful_math.R")
	output:
		out_combined = os.path.join(RESULTS_DIR, "plots", "enrichmentRecallCurves", "{comparison_name}.combinedCurves.pdf"),
		out_singles = os.path.join(RESULTS_DIR, "plots", "enrichmentRecallCurves", "{comparison_name}.individualCurves.pdf"),
		out_table = os.path.join(RESULTS_DIR, "plots", "enrichmentRecallCurves", "{comparison_name}.values.tsv.gz")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(ENV_DIR, "GWAS_env.yml")
	script:
		os.path.join(SCRIPTS_DIR, "visualization", "plot_enrichment_recall_curves.R")

# plot enrichment x recall and precision x recall scatterplots per comparison
rule plot_thresholded_performance_comparison:
	input: 
		variantOverlap = expand(os.path.join(RESULTS_DIR, "{method}", "variant_overlap", "enrichmentRecall.thresholded.traitByBiosample.tsv.gz"), method=config["methods"]),
		geneLinking = expand(os.path.join(RESULTS_DIR, "{method}", "gene_linking", "precisionRecall.traitByBiosample.tsv.gz"), method=config["methods"]),
		geneLinking_baseline = os.path.join(RESULTS_DIR, "baseline", "gene_linking", "precisionRecall.byTrait.tsv.gz"),
		traitGroups = os.path.join(RESULTS_DIR, "reference_configs", "trait_groups.tsv"),
		biosampleGroups =  os.path.join(RESULTS_DIR, "reference_configs", "biosample_groups.tsv"),
		allVariants = (os.path.join(RESULTS_DIR, "variants", "filteredGWASVariants.merged.sorted.tsv.gz")),
		colors = os.path.join(RESULTS_DIR, "plots", "colorPalette.tsv"),	
	params:
		biosample_names = lambda wildcards: comparisons.filter(pl.col("name")==wildcards.comparison_name)["biosample"].to_list(),
		trait_names = lambda wildcards: comparisons.filter(pl.col("name")==wildcards.comparison_name)["trait"].to_list(),
		p_threshold = config["thresholdPval"],
		helpful_math = os.path.join(SCRIPTS_DIR, "helpful_math.R")
	output:
		plots = os.path.join(RESULTS_DIR, "plots", "thresholdedPerformanceComparison", "{comparison_name}.scatter.pdf"),
		table_overlap = os.path.join(RESULTS_DIR, "plots", "thresholdedPerformanceComparison", "{comparison_name}.variantOverlap.tsv"),
		table_linking = os.path.join(RESULTS_DIR, "plots", "thresholdedPerformanceComparison", "{comparison_name}.geneLinking.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(ENV_DIR, "GWAS_env.yml")
	script:
		os.path.join(SCRIPTS_DIR, "visualization", "plot_thresholded_performance_comparison.R")