import os
import polars as pl

## PROCESS VARIANTS

# filter partition to distal noncoding (ABC, AllPeaks, Other, OtherIntron); filter background variants to distal noncoding
rule filter_sort_background_variants:
	input:
		partition = config["partition"],
		bgVariants = config["bgVariants"]
	params:
		chrSizes = config["chrSizes"]
	conda:
		os.path.join(ENV_DIR, "GWAS_env.yml")
	resources:
		mem_mb = determine_mem_mb
	output:
		partitionDistalNoncoding = temp(os.path.join(SCRATCH_DIR, "variants", "Partition.distalNoncoding.bed")),
		bgVarDistalNoncoding = temp(os.path.join(SCRATCH_DIR, "variants", "bgVariants.distalNoncoding.bed.gz"))
	shell:
		"""
		set +o pipefail;
		
		# filter partition to distal noncoding
		awk '$4=="ABC" || $4=="AllPeaks" || $4=="Other" || $4=="OtherIntron"' {input.partition} | bedtools sort -i stdin -faidx {params.chrSizes} > {output.partitionDistalNoncoding}

		# filter bg variants to distal noncoding
		cat {input.bgVariants} | bedtools sort -i stdin -faidx {params.chrSizes} | bedtools intersect -wa -sorted -a stdin -b {output.partitionDistalNoncoding} -g {params.chrSizes} | gzip > {output.bgVarDistalNoncoding}
		"""

# get count of distal noncoding background variants
rule background_variant_count:
	input:
		bgVarDistalNoncoding = os.path.join(SCRATCH_DIR, "variants", "bgVariants.distalNoncoding.bed.gz")
	output:
		bgVarCount = os.path.join(RESULTS_DIR, "variants", "bgVariants.distalNoncoding.count.txt")
	shell:
		"""
		set +o pipefail;
		zcat {input.bgVarDistalNoncoding} | wc -l > {output.bgVarCount}
		"""

# filter variant list to distal noncoding + by PIP and sort
# resultant columns:  chr,start,end,rsid,pip,CredibleSet,trait (with header; CredibleSet = "chr1:1000-2000")
rule filter_GWAS_variants:
	input: 
		variants = lambda wildcards: get_col2_from_col1(variant_key, "trait", wildcards.trait, "variant_file"),
		partitionDistalNoncoding = os.path.join(SCRATCH_DIR, "variants", "Partition.distalNoncoding.bed"),
	params:
		chrSizes = config["chrSizes"],
		thresholdPIP = config["thresholdPIP"]
	conda:
		os.path.join(ENV_DIR, "GWAS_env.yml")
	resources:
		mem_mb = determine_mem_mb
	output:
		variantsFiltered = temp(os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants", "{trait}.variantList.tsv"))
	shell:
		"""
		set +o pipefail;
		# header
		echo -e "chr\tstart\tend\trsid\tpip\tCredibleSet\ttrait" > {output.variantsFiltered}

		# rest of file
		cat {input.variants} | csvtk cut -t -f chr,start,end,rsid,pip,CredibleSet,Disease | sed 1d | awk '$5>{params.thresholdPIP}'  | \
			bedtools sort -i stdin -faidx {params.chrSizes} | bedtools intersect -wa -sorted -a stdin -b {input.partitionDistalNoncoding} -g {params.chrSizes} >> {output.variantsFiltered}
		"""

# concatenate all (filtered) variants into a single file and sort
# resultant columns:  chr,start,end,rsid,pip,CredibleSet,trait (with header; CredibleSet = "chr1:1000-2000")
rule combine_GWAS_variants:
	input:
		allVariants = expand(os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants", "{trait}.variantList.tsv"), trait=variant_key["trait"])
	params:
		variantDir = os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants"),
		chrSizes = config["chrSizes"],
	conda:
		os.path.join(ENV_DIR, "GWAS_env.yml")
	resources:
		mem_mb = determine_mem_mb
	output:
		header = temp(os.path.join(SCRATCH_DIR, "variants", "header.tsv")),
		variantsMerged = temp(os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants.merged.tsv")),
		variantsMergedSorted = temp(os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants.merged.sorted.tsv")),
		variantsMergedSortedGz = temp(os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants.merged.sorted.tsv.gz"))
	shell:
		"""
		set +o pipefail;

		directory={params.variantDir}
		first_file=$(find "$directory" -type f | head -n 1)
		head -n 1 "$first_file" > {output.header}

		# loop through all files in the directory
		for file in "$directory"/*; do
			tail -n +2 "$file" >> {output.variantsMerged}
		done

		# sort variants
		bedtools sort -i {output.variantsMerged} -faidx {params.chrSizes} > {output.variantsMergedSorted}

		cat {output.header} {output.variantsMergedSorted} | gzip > {output.variantsMergedSortedGz}

		"""