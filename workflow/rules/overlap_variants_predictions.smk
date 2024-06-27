import os
import polars as pl

## OVERLAP

# intersect background variants with thresholded predictions (trait by biosample) and save count (also for groups)
# resultant columns: biosample, bp_enhancers, bg_variant_count (with header)
rule compute_background_variant_overlap_and_enh_set_size:
	input:
		bgVarDistalNoncoding = os.path.join(SCRATCH_DIR, "variants", "bgVariants.distalNoncoding.bed.gz"),
		predictionsSorted = lambda wildcards: expand(os.path.join(SCRATCH_DIR, wildcards.method, "biosamples", "{biosample}", "enhancerPredictions.thresholded.bed.gz"), biosample=get_col2_from_col1(methods_config, "method", wildcards.method, "biosamples")),
		predictionsGroupsSorted = lambda wildcards: expand(os.path.join(SCRATCH_DIR, wildcards.method, "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.thresholded.merged.bed.gz"), biosampleGroup=get_col2_from_col1(methods_config, "method", wildcards.method, "biosampleGroups"))
	params:
		biosamples = lambda wildcards: get_col2_from_col1(methods_config, "method", wildcards.method, "biosamples"),
		biosampleGroups =  lambda wildcards: get_col2_from_col1(methods_config, "method", wildcards.method, "biosampleGroups"),
		chrSizes = config["chrSizes"]
	output:
		thresholdedPredictionStats = os.path.join(RESULTS_DIR, "{method}", "variant_overlap", "bgOverlap.enhancerSetSize.thresholdedPredictions.tsv")
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(ENV_DIR, "GWAS_env.yml")
	shell:
		"""			
		set +o pipefail;

		# create an output file with a header
		echo -e "biosample\tbp_enhancers\tbg_variant_count" > {output.thresholdedPredictionStats}

		# make array of biosamples + files
		IFS=' ' read -r -a biosampleArray <<< "{params.biosamples}"
		IFS=' ' read -r -a predFileArray <<< "{input.predictionsSorted}"

		# loop over the prediction files and corresponding biosamples
		for i in "${{!predFileArray[@]}}"
		do
			biosample=${{biosampleArray[$i]}}
			pred=${{predFileArray[$i]}}

			# calculate enh set size
			bps=$(zcat $pred | cut -f 1-3 | bedtools merge -i stdin | awk 'BEGIN {{FS=OFS="\t"}} {{print $3-$2}}' | awk '{{s+=$1}} END {{print s}}')
			# calculate bg variant overlap
			var_count=$(bedtools intersect -wa -sorted -g {params.chrSizes} -a <(zcat {input.bgVarDistalNoncoding}) -b <(zcat $pred)  | uniq | wc -l)

			# add to output file
			echo -e "${{biosample}}\t${{bps}}\t${{var_count}}" >> {output.thresholdedPredictionStats}
		done

		# repeat for groups
		IFS=' ' read -r -a biosampleGroupArray <<< "{params.biosampleGroups}"
		IFS=' ' read -r -a predGroupFileArray <<< "{input.predictionsGroupsSorted}"

		# loop over the prediction files and corresponding biosamples
		for i in "${{!predGroupFileArray[@]}}"
		do
			biosampleGroup=${{biosampleGroupArray[$i]}}
			predGroup=${{predGroupFileArray[$i]}}

			# calculate enh set size
			bps=$(zcat $predGroup | cut -f 1-3 | bedtools merge -i stdin | awk 'BEGIN {{FS=OFS="\t"}} {{print $3-$2}}' | awk '{{s+=$1}} END {{print s}}')
			# calculate bg variant overlap
			var_count=$(bedtools intersect -wa -sorted -g {params.chrSizes} -a <(zcat {input.bgVarDistalNoncoding}) -b <(zcat $predGroup)  | uniq | wc -l)

			# add to output file
			echo -e "${{biosampleGroup}}\t${{bps}}\t${{var_count}}" >> {output.thresholdedPredictionStats}
		done
		"""


# intersect variants from all traits with predictions (thresholded) from each biosample
# resultant columns with header: (variant 1-7) chr,start,end,rsid,pip,CredibleSet,trait; (prediction 8-13) chr,start,end,biosample,TargetGene,score
rule intersect_variants_predictions:
	input:
		variantsMerged =  os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants.merged.sorted.tsv.gz"),
		predictionsThresholded = os.path.join(SCRATCH_DIR, "{method}", "biosamples", "{biosample}", "enhancerPredictions.thresholded.bed.gz")
	params:
		chrSizes = config["chrSizes"]
	output:
		header =  temp(os.path.join(SCRATCH_DIR, "{method}", "biosamples", "{biosample}", "header.tsv")),
		variantsInt = temp(os.path.join(SCRATCH_DIR, "{method}", "biosamples", "{biosample}", "enhancerPredictions.variantIntersection.tsv")),
		variantsIntGz = os.path.join(SCRATCH_DIR, "{method}", "biosamples", "{biosample}", "enhancerPredictions.variantIntersection.tsv.gz")
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(ENV_DIR, "GWAS_env.yml")
	shell:
		"""
		set +o pipefail;

		echo -e "varChr\tvarStart\tvarEnd\trsid\tPIP\tCredibleSet\ttrait\tpredChr\tpredStart\tpredEnd\tbiosample\tTargetGene\tpredScore" > {output.header}
		bedtools intersect -wa -wb -sorted -g {params.chrSizes} -a <(zcat {input.variantsMerged}) -b <(zcat {input.predictionsThresholded}) > {output.variantsInt}

		cat {output.header} {output.variantsInt} | gzip > {output.variantsIntGz}
		"""

# intersect merged predictions from biosample groups with variants
# resultant columns with header: (variant 1-7) chr,start,end,rsid,pip,CredibleSet,trait; (prediction 8-10) chr,start,end
rule intersect_variants_merged_predictions:
	input:
		variantsMerged = os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants.merged.sorted.tsv.gz"),
		predictionsThresholded = os.path.join(RESULTS_DIR, "{method}", "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.thresholded.merged.bed.gz")
	params:
		chrSizes = config["chrSizes"]
	output:
		header =  temp(os.path.join(SCRATCH_DIR, "{method}", "biosamples", "{biosampleGroup}", "header.tsv")),
		variantsInt = temp(os.path.join(SCRATCH_DIR, "{method}", "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.variantIntersection.tsv")),
		variantsIntGz = os.path.join(SCRATCH_DIR, "{method}", "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.variantIntersection.tsv.gz")
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(ENV_DIR, "GWAS_env.yml")
	shell:
		"""
		set +o pipefail;

		echo -e "varChr\tvarStart\tvarEnd\trsid\tPIP\tCredibleSet\ttrait\tpredChr\tpredStart\tpredEnd" > {output.header}
		bedtools intersect -wa -wb -sorted -g {params.chrSizes} -a <(zcat {input.variantsMerged}) -b <(zcat {input.predictionsThresholded}) > {output.variantsInt}

		cat {output.header} {output.variantsInt} | gzip > {output.variantsIntGz}
		"""

rule calculate_enrichment_recall:
	input:
		variantsInt = lambda wildcards: expand(os.path.join(SCRATCH_DIR, wildcards.method, "biosamples", "{biosample}", "enhancerPredictions.variantIntersection.tsv.gz"), biosample=get_col2_from_col1(methods_config, "method", wildcards.method, "biosamples")),
		variantsIntGroups = lambda wildcards: expand(os.path.join(SCRATCH_DIR, wildcards.method, "biosampleGroups", "{biosampleGroup}", "enhancerPredictions.variantIntersection.tsv.gz"), biosampleGroup=get_col2_from_col1(methods_config, "method", wildcards.method, "biosampleGroups")),
		thresholdedPredictionStats = os.path.join(RESULTS_DIR, "{method}", "variant_overlap", "bgOverlap.enhancerSetSize.thresholdedPredictions.tsv"),
		variantsMerged = os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants.merged.sorted.tsv.gz"),
		bgVarCount = os.path.join(RESULTS_DIR, "variants", "bgVariants.distalNoncoding.count.txt")
	params:
		biosamples = lambda wildcards: get_col2_from_col1(methods_config, "method", wildcards.method, "biosamples"),
		biosampleGroups = lambda wildcards: get_col2_from_col1(methods_config, "method", wildcards.method, "biosampleGroups"),
		thresholdPval = config["thresholdPval"]
	output:
		res = os.path.join(RESULTS_DIR, "{method}", "variant_overlap", "enrichmentRecall.traitByBiosample.tsv.gz")
	resources:
		mem_mb = determine_mem_mb
	conda: 
		os.path.join(ENV_DIR, "GWAS_env.yml")
	script:
		os.path.join(SCRIPTS_DIR, "enrichment_recall", "calculate_enrichment_recall.R")


