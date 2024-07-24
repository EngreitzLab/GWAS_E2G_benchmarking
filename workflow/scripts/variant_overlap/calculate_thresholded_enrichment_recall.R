suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  source(snakemake@params$helpful_math)})

## files from snakemake
varInt_files = snakemake@input$variantsInt %>% strsplit(" ") %>% unlist()
varIntGroup_files = snakemake@input$variantsIntGroups %>% strsplit(" ") %>% unlist()
biosamples = snakemake@params$biosamples %>% strsplit(" ") %>% unlist()
biosampleGroups = snakemake@params$biosampleGroups %>% strsplit(" ") %>% unlist()
traitGroups = snakemake@params$traitGroups %>% strsplit(" ") %>% unlist()
print(traitGroups)
predStats_file = snakemake@input$thresholdedPredictionStats
GWASVar_file = snakemake@input$variantsMerged
bgVarCount = readLines(snakemake@input$bgVarCount) %>% as.numeric()
p_threshold = snakemake@params$thresholdPval %>% as.numeric()
method_this = snakemake@wildcards$method
out_file = snakemake@output$res

## create dataframe with the following columns: biosample[/group],trait,nVariantsOverlappingEnhancers,nVariantsTotal,nCommonVariantsOverlappingEnhancers,nCommonVariantsTotal,group(bool)

# get nVariantsTotal  (number of variants per trait)
var = fread(GWASVar_file, sep="\t", header=TRUE) # chr,start,end,rsid,CredibleSet,trait 
var = dplyr::select(var, -c(rsid, CredibleSet)) %>% distinct() %>% group_by(trait) %>% summarize(nVariantsTotal=n())
print(var)

# get nCommonVariantsoverlappingEnhancers
bgVar = fread(predStats_file, sep="\t", header=TRUE) %>% setNames(c("biosample", "bpEnhancers", "nCommonVariantsOverlappingEnhancers"))

# get nVariantsOverlappingEnhancers
varInt_files_all = c(varInt_files, varIntGroup_files)
biosamples_all = c(biosamples, biosampleGroups)

for (i in 1:length(biosamples_all)){
	biosample_this = biosamples_all[i]
	varInt_this = fread(varInt_files_all[i], sep="\t", header=TRUE)  #(variant 1-6) chr,start,end,rsid,CredibleSet,trait; (prediction 7-12) chr,start,end,biosample,TargetGene,score
	varInt_counts = dplyr::select(varInt_this, varChr,varStart,varEnd,trait) %>% distinct() %>%
		group_by(trait) %>% summarize(nVariantsOverlappingEnhancers=n()) %>% 
		mutate(biosample = biosamples_all[i])

	if (i==1) {df = varInt_counts} else {df = rbind(df, varInt_counts)}
}

# assemble
df = left_join(df, var, by=c("trait")) %>% # add nVariantsTotal
	left_join(bgVar, by=c("biosample")) %>% # add bpEnhancers, nCommonVariantsOverlappingEnhancers
	mutate(nCommonVariantsTotal = bgVarCount,
		recall = nVariantsOverlappingEnhancers/nVariantsTotal,
		enrichment = recall/(nCommonVariantsOverlappingEnhancers/nCommonVariantsTotal),
		method = method_this,
		biosampleGroup = ifelse(biosample %in% biosampleGroups, TRUE, FALSE),
		traitGroup = ifelse(trait %in% traitGroups, TRUE, FALSE)) 
df = add_enrichment_statistics(df, p_threshold)
df = add_recall_overlap_statistics(df, p_threshold)

fwrite(df, file=out_file, sep="\t", quote=F, row.names=F, col.names=T)
