import os
import polars as pl

## LINKING VARIANTS TO GENES

# intersect variants from all traits with predictions (thresholded) from each biosample
# resultant columns with header: (variant 1-7) chr,start,end,rsid,pip,CredibleSet,trait; (prediction 8-13) chr,start,end,biosample,TargetGene,score
rule intersect_variants_concat_predictions_with_gene:
	input:
		variantsMerged = ancient(os.path.join(RESULTS_DIR, "variants", "filteredGWASVariants.merged.sorted.tsv.gz")),
		predictionsConcatThresholded = os.path.join(SCRATCH_DIR, "{method}", "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.thresholded.concat.bed.gz")
	params:
		chrSizes = config["chrSizes"]
	output:
		header =  temp(os.path.join(SCRATCH_DIR, "{method}", "biosamples", "{biosampleGroup}", "headerWithGene.tsv")),
		variantsInt = temp(os.path.join(SCRATCH_DIR, "{method}", "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.variantIntersection.withGene.tsv")),
		variantsIntGz = os.path.join(SCRATCH_DIR, "{method}", "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.variantIntersection.withGene.tsv.gz")
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(ENV_DIR, "GWAS_env.yml")
	shell:
		"""
		set +o pipefail;

		echo -e "varChr\tvarStart\tvarEnd\trsid\tCredibleSet\ttrait\tpredChr\tpredStart\tpredEnd\tbiosample\tTargetGene\tpredScore" > {output.header}
		bedtools intersect -wa -wb -sorted -g {params.chrSizes} -a <(zcat {input.variantsMerged} | sed 1d) -b <(zcat {input.predictionsConcatThresholded}) > {output.variantsInt}

		cat {output.header} {output.variantsInt} | gzip > {output.variantsIntGz}
		"""

# calculate precision and recall of linking variants to their target gene  
rule evalute_gene_linking:
	input:
		variantsInt = lambda wildcards: expand(os.path.join(SCRATCH_DIR, wildcards.method, "biosamples", "{biosample}", "enhancerPredictions.variantIntersection.tsv.gz"), biosample=get_col2_from_col1(methods_config, "method", wildcards.method, "biosamples")),
		variantsIntGroups = lambda wildcards: expand(os.path.join(SCRATCH_DIR, wildcards.method, "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.variantIntersection.withGene.tsv.gz"), biosampleGroup=get_col2_from_col1(methods_config, "method", wildcards.method, "biosampleGroups")),
	params:
		TSS = config["TSS"],
		biosamples = lambda wildcards: get_col2_from_col1(methods_config, "method", wildcards.method, "biosamples"),
		biosampleGroups = lambda wildcards: get_col2_from_col1(methods_config, "method", wildcards.method, "biosampleGroups"),
		genePrioritizationTable = config["genePrioritizationTable"], 
		thresholdPval = config["thresholdPval"],
		numPoPSGenes = config["numPoPSGenes"],
		numPredGenes = config["numPredGenes"],
		boolean = lambda wildcards: get_col2_from_col1(methods_config, "method", wildcards.method, "boolean"),
		helpful_math = os.path.join(SCRIPTS_DIR, "helpful_math.R")
	output:
		res = os.path.join(RESULTS_DIR, "{method}", "gene_linking", "precisionRecall.traitByBiosample.tsv.gz")
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(ENV_DIR, "GWAS_env.yml")
	script:
		os.path.join(SCRIPTS_DIR, "gene_linking", "evaluate_gene_linking.R")

rule evaluate_baseline_predictors_gene_linking:
	params:
		TSS = config["TSS"],
		genePrioritizationTable = config["genePrioritizationTable"],
		numPoPSGenes = config["numPoPSGenes"],
		thresholdPval = config["thresholdPval"],
		helpful_math = os.path.join(SCRIPTS_DIR, "helpful_math.R")
	output:
		res = os.path.join(RESULTS_DIR, "baseline", "gene_linking", "precisionRecall.byTrait.tsv.gz")
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(ENV_DIR, "GWAS_env.yml")
	script:
		os.path.join(SCRIPTS_DIR, "gene_linking", "evaluate_baseline_predictors.R")
