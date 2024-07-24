import os
import polars as pl

## PROCESS PREDICTIONS

# initial processing: read in, sort, adjust score, filter to gene universe
# resultant columns:  chr,start,end,biosample,TargetGene,score (no header)
rule process_predictions:
	input:
		predFile = lambda wildcards: get_pred_file(wildcards.method, wildcards.biosample),
		geneUniverse = config["genes"]
	params:
		scriptsDir = SCRIPTS_DIR,
		chrSizes = config["chrSizes"],
		scoreCol = lambda wildcards: get_col2_from_col1(methods_config, "method", wildcards.method, "score_col"),
		inversePred = lambda wildcards: get_col2_from_col1(methods_config, "method", wildcards.method, "inverse_predictor"),
	output:
		predictionsTemp = temp(os.path.join(SCRATCH_DIR, "{method}", "biosamples", "{biosample}", "enhancerPredictions.temp.tsv")),
		predictionsSorted = (os.path.join(SCRATCH_DIR, "{method}", "biosamples", "{biosample}", "enhancerPredictions.sorted.bed.gz"))
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(ENV_DIR, "GWAS_env.yml")
	shell:
		"""			
		# sort predictions file: remove # from header,select columns,remove header, remove rows with blanks
		set +o pipefail;
		if [[ {input.predFile} == *.gz ]]
		then
			zcat {input.predFile} | awk 'NR==1{{sub(/^#*/, "")}}1' | csvtk cut -t -f chr,start,end,TargetGene,{params.scoreCol} | sed 1d | awk 'NF==5{{print}}{{}}' | bedtools sort -i stdin -faidx {params.chrSizes} > {output.predictionsTemp}
		else
			cat {input.predFile} | awk 'NR==1{{sub(/^#*/, "")}}1' | csvtk cut -t -f chr,start,end,TargetGene,{params.scoreCol} | sed 1d | awk 'NF==5{{print}}{{}}' | bedtools sort -i stdin -faidx {params.chrSizes} > {output.predictionsTemp}
		fi

		# invert score if inverted predictor and filter to gene universe and set biosample column as 4th 
		Rscript {params.scriptsDir}/preprocessing/process_predictions.R --input {output.predictionsTemp}  --genes {input.geneUniverse} --biosample {wildcards.biosample} --invert {params.inversePred}  | gzip > {output.predictionsSorted}
			
		"""

# threshold predictions (columns: chr,start,end,biosample,TargetGene,score; no header)
rule threshold_predictions: 
	input:
		predictionsSorted = os.path.join(SCRATCH_DIR, "{method}", "biosamples", "{biosample}", "enhancerPredictions.sorted.bed.gz")
	params:
		threshold = lambda wildcards: get_col2_from_col1(methods_config, "method", wildcards.method, "threshold")
	conda: 
		os.path.join(ENV_DIR, "GWAS_env.yml")
	resources:
		mem_mb = determine_mem_mb
	output:
		predictionsThresholded = os.path.join(SCRATCH_DIR, "{method}", "biosamples", "{biosample}", "enhancerPredictions.thresholded.bed.gz")
	shell:
		"""
		set +o pipefail;

		zcat {input.predictionsSorted} | awk '$6>={params.threshold}' | gzip > {output.predictionsThresholded}
		"""

# concatenate thresholded predictions for a biosample group (no header, columns: chr,start,end,biosample,TargetGene,score)
rule concatenate_thresholded_predictions: 
	input:
		predictionsThresholded = lambda wildcards: expand(os.path.join(SCRATCH_DIR, wildcards.method, "biosamples", "{biosample}", "enhancerPredictions.thresholded.bed.gz"), biosample=get_single_biosamples(wildcards.method, wildcards.biosampleGroup))
	params:
		chrSizes = config["chrSizes"]
	conda: 
		os.path.join(ENV_DIR, "GWAS_env.yml")
	resources:
		mem_mb = determine_mem_mb
	output:
		predictionsConcat = temp(os.path.join(SCRATCH_DIR, "{method}", "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.thresholded.concat.bed")),
		predictionsConcatGz = os.path.join(SCRATCH_DIR, "{method}", "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.thresholded.concat.bed.gz")
	shell: 
			"""			
			set +o pipefail;

			# make array of files
			IFS=' ' read -r -a predFileArray <<< "{input.predictionsThresholded}"

			# loop over the prediction files
			for i in "${{!predFileArray[@]}}"
			do
				pred=${{predFileArray[$i]}}
				zcat $pred >> {output.predictionsConcat}
			done

			cat {output.predictionsConcat} | bedtools sort -i stdin -faidx {params.chrSizes} | gzip > {output.predictionsConcatGz}
			"""

# merge regions for biosample group (columns: chr,start,end with no header)
rule merge_thresholded_predictions:
	input:
		predictionsConcat = os.path.join(SCRATCH_DIR, "{method}", "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.thresholded.concat.bed.gz")
	params:
		chrSizes = config["chrSizes"]
	conda: 
		os.path.join(ENV_DIR, "GWAS_env.yml")
	resources:
		mem_mb = determine_mem_mb
	output:
		predictionsMerged = os.path.join(SCRATCH_DIR, "{method}", "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.thresholded.merged.bed.gz")
	shell:
		"""
		set +o pipefail;

		zcat {input.predictionsConcat} | cut -f 1-3 | bedtools sort -i stdin -faidx {params.chrSizes} | bedtools merge -i stdin | gzip > {output.predictionsMerged}
		"""

# generate threshold span based on quantiles of interescting predictions (per method)
rule generate_quantile_threshold_span:
	input:
		varInt = lambda wildcards: [os.path.join(SCRATCH_DIR, wildcards.method, "biosamples", biosample, "enhancerPredictions.variantIntersection.tsv.gz") for biosample in get_comparison_biosamples_per_method(comparisons, wildcards.method)[0]]
	params:
		nSteps = config["nThresholdSteps"],
		comparisonsTable = config["comparisonsTable"],
		threshold = lambda wildcards: get_col2_from_col1(methods_config, "method", wildcards.method, "threshold"),
		biosamples = lambda wildcards: get_comparison_biosamples_per_method(comparisons, wildcards.method)[1],
		traits = lambda wildcards: get_comparison_biosamples_per_method(comparisons, wildcards.method)[2],
		boolean = lambda wildcards: get_col2_from_col1(methods_config, "method", wildcards.method, "boolean")
	output:
		thresholdSpan = os.path.join(RESULTS_DIR, "{method}", "thresholdSpan.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda:
		os.path.join(ENV_DIR, "GWAS_env.yml")
	script:
		os.path.join(SCRIPTS_DIR, "preprocessing", "generate_quantile_threshold_span.R")
		