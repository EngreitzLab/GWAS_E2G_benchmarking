suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
})


make_bed12_row <- function(df, CS_id){
	df_this = dplyr::filter(df, CS_trait==CS_id) %>%
		arrange(start) %>%
		mutate(rel_start=start-min(start)) # this might be off by one
	
	nSNP = nrow(df_this)
	trait = df_this$trait[1]

	res = with(df_this, c(chr=chr[1], start=min(start), end=max(end), name=CS_trait[1],
		score=nSNP, strand=".", thickStart=min(start), thickEnd=max(end),  itemRgb="0,0,0",
		blockCount=nSNP, blockSizes=paste0(rep("1,",nSNP), collapse=""),
		blockStarts=paste0(rel_start,collapse=",")))

	res = data.frame(as.list(res))
	return(res)
}

process_merged_variants <- function(input_file, traits_to_remove, output_file){
	# header: chr	start	end	rsid	CredibleSet	trait
	# CredibleSet format: chr1:0-4889757-2; might be in hg19

	df = fread(input_file, header=TRUE, sep="\t") %>%
		dplyr::filter(!(trait %in% traits_to_remove)) %>%
		mutate(CS_trait  = paste(trait, CredibleSet, sep="_")) 

	CS_vector = unique(df$CS_trait)
	nCS = length(CS_vector)
	bed12_list = vector("list", nCS)

	for (i in 1:nCS){
		bed12_list[[i]] = make_bed12_row(df, CS_vector[i])
	}

	bed12 = rbindlist(bed12_list) %>% as_tibble()
	print(head(bed12))

	fwrite(bed12, output_file, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, scipen=999)

}






### MAIN
input_file = "/oak/stanford/groups/engreitz/Users/sheth/GWAS_benchmarking_working/GWAS_E2G_benchmarking/results/2024_0711_sc_benchmarking_v1/variants/filteredGWASVariants.merged.sorted.tsv.gz"
output_file = "/oak/stanford/groups/engreitz/igv/hg38/scE2G/UKBB_GWAS/merged_variants.bed"
traits_to_remove = c("ALL", "RBC_traits", "WBC_traits", "Metabolic_traits", "Autoimmune_traits", "Blood_counts", "T2D_related")

# run
process_merged_variants(input_file, traits_to_remove, output_file)