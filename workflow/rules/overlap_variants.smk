import os
import pandas as pd

# function returning list of files from checkpoint
def get_chrom_files(wildcards):
	with checkpoints.define_chromosomes.get(method=wildcards.method, tissue=wildcards.tissue).output.chromsomesPresent.open() as file:
		bim_chr = [f"chr{num}" for num in range(1, 23)] # chr1, chr2, ..., chr22
		chr_list = [line.rstrip() for line in file if line.rstrip() in bim_chr]
		file_list = [os.path.join(config["outDir"], method, "tissues", tissue, "bgVariantAnnotations", f"{chrom}.bed.gz") for chrom in chr_list]
		print(file_list)
		return file_list

def get_variant_files(trait_group, what_to_return):
	file_list = []
	name_list = []
	for index, row in variant_key.iterrows():
		groups_this = [x.strip() for x in row["trait_group"].split(',')] 
		if trait_group in groups_this:
			file_list.append(row["trait_file"])
			name_list.append(row["trait"])
	if what_to_return=="files": 
		return file_list
	else:
		return name_list

rule overlap_predictions_with_variants: # for GM12878 (WBC) and K562 (RBC)
	input:
		bgVariantAnnotations = get_chrom_files # based on chromosomes from checkpoint
	params:
		finemappedVarFiles = lambda wildcards: get_variant_files(wildcards.trait_group, "files"),
		finemappedVarTraits = lambda wildcards: get_variant_files(wildcards.trait_group, "traits"),
		nBoostrap = 100
	conda: 
		os.path.join(config["envDir"], "GWAS_env.yml")
	resources:
		mem_mb = determine_mem_mb
	output:
		results = os.path.join(config["outDir"], "{method}", "tissues", "{tissue}", "{trait_group}_overlap_results.tsv")
	script:
		os.path.join(config["codeDir"], "overlap_finemapped_variants.R")