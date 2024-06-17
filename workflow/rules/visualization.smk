import os
import pandas as pd


rule plot_performance: # per comparison
	input:
		overlaps = [os.path.join(config["outDir"], method, "tissues", "{tissue}", "{trait_group}_overlap_results.tsv") for method in config["methods"]]
	params:
		methods = config["methods"],
		methodsConfig = config["methodsTable"] # to get colors
	conda: 
		os.path.join(config["envDir"], "GWAS_env.yml")
	resources:
		mem_mb = determine_mem_mb
	output:
		outPlot = os.path.join(config["outDir"], "plots",  "Tissue{tissue}_Trait{trait_group}_overlap_results.pdf"),
		outFile = os.path.join(config["outDir"], "plots",  "Tissue{tissue}_Trait{trait_group}_overlap_results.tsv")
	script:
		os.path.join(config["codeDir"], "plot_performance_metrics.R")