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
# resultant columns:  chr,start,end,rsid,CredibleSet,trait (with header; CredibleSet = "chr1:1000-2000")
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
		variantsFiltered = os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants", "{trait}.variantList.tsv")
	shell:
		"""
		set +o pipefail;
		# header
		echo -e "chr\tstart\tend\trsid\tCredibleSet\ttrait" > {output.variantsFiltered}

		# rest of file
		cat {input.variants} | csvtk cut -t -f chr,start,end,rsid,pip,CredibleSet,trait | sed 1d | awk '$5>{params.thresholdPIP}'  | cut -f1,2,3,4,6,7 | \
			bedtools sort -i stdin -faidx {params.chrSizes} | bedtools intersect -wa -sorted -a stdin -b {input.partitionDistalNoncoding} -g {params.chrSizes} >> {output.variantsFiltered}
		"""

# merge & deduplicate variant lists for trait groups
rule merge_GWAS_variants_trait_group:
	input: 
		variantFiles = lambda wildcards: expand((os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants", "{trait}.variantList.tsv")), trait=config["traitGroups"][wildcards.traitGroup])
	params:
		chrSizes = config["chrSizes"],
	conda:
		os.path.join(ENV_DIR, "GWAS_env.yml")
	resources:
		mem_mb = determine_mem_mb
	output:
		variantsTraitGroupConcat =temp(os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants", "traitGroups", "{traitGroup}.variantList.concat.tsv")),
		variantsTraitGroup = os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants", "traitGroups", "{traitGroup}.variantList.tsv")
	shell:
		"""
		set +o pipefail;
		# header
		echo -e "chr\tstart\tend\trsid\tCredibleSet\ttrait" >> {output.variantsTraitGroup}

		# make array of files
		IFS=' ' read -r -a variantFileArray <<< "{input.variantFiles}"

		# loop over the variant files and concat 
		for i in "${{!variantFileArray[@]}}"
		do
			varFile=${{variantFileArray[$i]}}
			# remove header and trait column, add to concat list
			cat $varFile | sed 1d | cut -f1-5 >> {output.variantsTraitGroupConcat}
		done

		# sort/dedup and add trait column, add to final file
		wc -l {output.variantsTraitGroupConcat}
		cat {output.variantsTraitGroupConcat} | sort | uniq | sed 's/$/\t{wildcards.traitGroup}/' >> {output.variantsTraitGroup}
		wc -l {output.variantsTraitGroup}
		"""
	

# concatenate all (filtered) variants into a single file and sort
# resultant columns:  chr,start,end,rsid,CredibleSet,trait (with header; CredibleSet = "chr1:1000-2000")
rule combine_GWAS_variants:
	input:
		allTraits = expand(os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants", "{trait}.variantList.tsv"), trait=variant_key["trait"]),
		allTraitGroups = expand(os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants", "traitGroups", "{traitGroup}.variantList.tsv"), traitGroup=[*config["traitGroups"]])
	params:
		chrSizes = config["chrSizes"],
	conda:
		os.path.join(ENV_DIR, "GWAS_env.yml")
	resources:
		mem_mb = determine_mem_mb
	output:
		header = temp(os.path.join(SCRATCH_DIR, "variants", "header.tsv")),
		variantsMerged = temp(os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants.merged.tsv")),
		variantsMergedSorted = temp(os.path.join(SCRATCH_DIR, "variants", "filteredGWASVariants.merged.sorted.tsv")),
		variantsMergedSortedGz = (os.path.join(RESULTS_DIR, "variants", "filteredGWASVariants.merged.sorted.tsv.gz"))
	shell:
		"""
		set +o pipefail;

		# header
		echo -e "chr\tstart\tend\trsid\tCredibleSet\ttrait" > {output.header}

		# arrays of files
		IFS=' ' read -r -a traitFileArray <<< "{input.allTraits}"
		IFS=' ' read -r -a traitGroupFileArray <<< "{input.allTraitGroups}"

		# loop over trait files and concat 
		for i in "${{!traitFileArray[@]}}"
		do
			traitFile=${{traitFileArray[$i]}}
			# remove header, add to concat list
			cat $traitFile | sed 1d >> {output.variantsMerged}
		done
		
		# add trait groups
		for i in "${{!traitGroupFileArray[@]}}"
		do
			traitGroupFile=${{traitGroupFileArray[$i]}}
			# remove header, add to concat list
			cat $traitGroupFile | sed 1d >> {output.variantsMerged}
		done

		# sort variants
		bedtools sort -i {output.variantsMerged} -faidx {params.chrSizes} > {output.variantsMergedSorted}

		cat {output.header} {output.variantsMergedSorted} | gzip > {output.variantsMergedSortedGz}

		"""