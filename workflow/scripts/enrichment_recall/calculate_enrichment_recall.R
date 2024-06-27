suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

## files from snakemake
varInt_files = snakemake@input$variantsInt %>% strsplit(" ") %>% unlist()
varIntGroup_files = snakemake@input$variantsIntGroups %>% strsplit(" ") %>% unlist()
biosamples = snakemake@params$biosamples %>% strsplit(" ") %>% unlist()
groups = snakemake@params$biosampleGroups %>% strsplit(" ") %>% unlist()
predStats_file = snakemake@input$thresholdedPredictionStats
GWASVar_file = snakemake@input$variantsMerged
bgVarCount = readLines(snakemake@input$bgVarCount) %>% as.numeric()
p_threshold = snakemake@params$thresholdPval %>% as.numeric()
method_this = snakemake@wildcards$method
out_file = snakemake@output$res

## create dataframe with the following columns: biosample[/group],trait,nVariantsOverlappingEnhancers,nVariantsTotal,nCommonVariantsOverlappingEnhancers,nCommonVariantsTotal,group(bool)
# get nVariantsTotal  (number of variants per trait)
var = fread(GWASVar_file, sep="\t", header=TRUE) # chr,start,end,rsid,pip,CredibleSet,trait 
var = dplyr::select(var, -c(rsid, pip, CredibleSet)) %>% distinct() %>% group_by(trait) %>% summarize(nVariantsTotal=n())

# get nCommonVariantsoverlappingEnhancers
bgVar = fread(predStats_file, sep="\t", header=TRUE) %>% setNames(c("biosample", "bpEnhancers", "nCommonVariantsOverlappingEnhancers"))

# get nVariantsOverlappingEnhancers
varInt_files_all = c(varInt_files, varIntGroup_files)
biosamples_all = c(biosamples, groups)

for (i in 1:length(biosamples_all)){
	biosample_this = biosamples_all[i]
	varInt_this = fread(varInt_files_all[i], sep="\t", header=TRUE)  #(variant 1-7) chr,start,end,rsid,pip,CredibleSet,trait; (prediction 8-13) chr,start,end,biosample,TargetGene,score
	varInt_this = dplyr::select(varInt_this, varChr,varStart,varEnd,trait) %>% distinct() %>%
		group_by(trait) %>% summarize(nVariantsOverlappingEnhancers=n()) %>% 
		mutate(biosample = biosamples_all[i], group=ifelse(i>length(biosamples), TRUE, FALSE))

	if (i==1) {df = varInt_this} else {df = rbind(df, varInt_this)}
}

# assemble
df = left_join(df, var, by=c("trait")) %>% # add nVariantsTotal
	left_join(bgVar, by=c("biosample")) %>% # add bpEnhancers, nCommonVariantsOverlappingEnhancers
	mutate(nCommonVariantsTotal = bgVarCount,
		recall = nVariantsOverlappingEnhancers/nVariantsTotal,
		enrichment = recall/(nCommonVariantsOverlappingEnhancers/nCommonVariantsTotal),
		method = method_this) 

## statistics on enrichment + recall
z = qnorm(p_threshold/2, lower.tail=FALSE) # e.g. 1.96

# CI + pval (hypergeometric) for enrichment aka risk ratio > 1
# CI ref: https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_confidence_intervals/bs704_confidence_intervals8.html
df = mutate(df, SE_log_enr  =sqrt(((nVariantsTotal-nVariantsOverlappingEnhancers)/nVariantsOverlappingEnhancers)/nVariantsTotal) + ((nCommonVariantsTotal-nCommonVariantsOverlappingEnhancers)/nCommonVariantsOverlappingEnhancers)/nCommonVariantsTotal,
	CI_enr_low = exp(log(enrichment)-z*SE_log_enr),
	CI_enr_high = exp(log(enrichment)+z*SE_log_enr),
	p_enr = phyper(nVariantsOverlappingEnhancers, nVariantsTotal, nCommonVariantsTotal, nVariantsOverlappingEnhancers + nCommonVariantsOverlappingEnhancers, log.p=FALSE, lower.tail=FALSE),
	p_adjust_enr = p.adjust(p_enr, method="bonferroni"))

# CI for recall (Wilson score interval): https://users.stat.ufl.edu/~aa/articles/agresti_coull_1998.pdf
df = mutate(df, recall_adjust = (nVariantsOverlappingEnhancers+2)/(nVariantsTotal+4),
	CI_recall_low = recall_adjust - z*sqrt(recall_adjust * (1-recall_adjust)/(nVariantsTotal+4)),
	CI_recall_high = recall_adjust + z*sqrt(recall_adjust * (1-recall_adjust)/(nVariantsTotal+4)))

fwrite(df, file=out_file, sep="\t", quote=F, row.names=F, col.names=T)
