suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  source(snakemake@params$helpful_math)})

### files from snakemake
varInt_files = snakemake@input$variantsInt %>% strsplit(" ")
bgVarInt_files = snakemake@input$bgVariantsInt %>% strsplit(" ")
GWASVar_file = snakemake@input$variantsMerged
thresholdSpan = fread(snakemake@input$thresholdSpan, sep="\t")
bgVarCount_total = readLines(snakemake@input$bgVarCount) %>% as.numeric()
biosamples = snakemake@params$biosamples %>% strsplit(" ") %>% unlist()
traitGroups = snakemake@params$traitGroups %>% strsplit(" ") %>% unlist()
p_threshold = snakemake@params$thresholdPval %>% as.numeric()
method_this = snakemake@wildcards$method
biosample_name = snakemake@wildcards$biosample # group or single name
out_file = snakemake@output$res 

### create dataframe with the following columns: threshold,biosample[/group],trait,nVariantsOverlappingEnhancers,nVariantsTotal,nCommonVariantsOverlappingEnhancers,nCommonVariantsTotal

## get nVariantsTotal  (number of variants per trait)
var = fread(GWASVar_file, sep="\t", header=TRUE) # chr,start,end,rsid,CredibleSet,trait 
var = dplyr::select(var, -c(rsid, CredibleSet)) %>% distinct() %>% group_by(trait) %>% summarize(nVariantsTotal=n())

## get nVariantsOverlappingEnhancers and nCommonVariantsOverlappingEnhancers
# read in + format files
varInt = lapply(varInt_files, fread, sep="\t") %>% #var[varChr, varStart, varEnd, rsid, CredibleSet, trait],  pred[predChr, predStart, predEnd, biosample, TargetGene, predScore]
	rbindlist() %>% as_tibble() %>%
	dplyr::select(varChr, varStart, varEnd, trait, predScore) # everything else is irrelevant?
bgVarInt = lapply(bgVarInt_files, fread, sep="\t") %>% # bgvar[varChr, varStart, varEnd, rsid], pred[predChr, predStart, predEnd, biosample, TargetGene, predScore]
	rbindlist() %>% as_tibble() %>%
	dplyr::select(varChr, varStart, varEnd, predScore)

res_list = vector("list", nrow(thresholdSpan))
for (i in 1:nrow(thresholdSpan)){
	threshold_this = thresholdSpan$threshold[i]

	bg_count = dplyr::filter(bgVarInt, predScore>=threshold_this) %>%
		dplyr::select(-predScore) %>% distinct() %>% nrow()

	res_list[[i]] = dplyr::filter(varInt, predScore>=threshold_this) %>% 
		dplyr::select(-predScore) %>% distinct() %>%
		group_by(trait) %>% summarize(nVariantsOverlappingEnhancers=n()) %>%
		mutate(nCommonVariantsOverlappingEnhancers=bg_count, threshold=threshold_this)
}

# assemble
res = rbindlist(res_list) %>% as_tibble() %>%
	left_join(var, by=c("trait")) %>% # add nVariantsTotal
	mutate(nCommonVariantsTotal = bgVarCount_total,
		recall = nVariantsOverlappingEnhancers/nVariantsTotal,
		enrichment = recall/(nCommonVariantsOverlappingEnhancers/nCommonVariantsTotal),
		biosample = biosample_name,
		method = method_this,
		biosampleGroup = ifelse(length(biosamples)>1, TRUE, FALSE),
		traitGroup = ifelse(trait %in% traitGroups, TRUE, FALSE))

res = add_enrichment_statistics(res, p_threshold)
res = add_recall_overlap_statistics(res, p_threshold)

fwrite(res, file=out_file, sep="\t", quote=F, row.names=F, col.names=T)
