suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  source(snakemake@params$helpful_math)})

## files from snakemake
genePri_file = snakemake@params$genePrioritizationTable
p_threshold = snakemake@params$thresholdPval %>% as.numeric()
numPoPSGenes = snakemake@params$numPoPSGenes %>% as.numeric()
out_file = snakemake@output$res

#### HELPER FUNCTIONS
# return a df row with results & correct column names
get_results_vector  <- function(df, trait_this, method){
	df = dplyr::filter(df, trait==trait_this) # CredibleSet, trait, TargetGene, TruthGene

	nCredibleSetsWithPredictions = dplyr::select(df, CredibleSet, TruthGene) %>% distinct() %>% nrow() # unique CS-gene pairs with predictions (should be all)
	nCredibleSetGenePredictions = dplyr::select(df, CredibleSet, TargetGene) %>% distinct() %>% nrow()
	nCredibleSetGenePredictionsCorrect = dplyr::filter(df, TargetGene==TruthGene) %>% distinct() %>% nrow()		

	res = c(trait_this, method, nCredibleSetsWithPredictions, nCredibleSetGenePredictions, nCredibleSetGenePredictionsCorrect)
	names(res) = c("trait", "method", "nCredibleSetsWithPredictions", "nCredibleSetGenePredictions", "nCredibleSetGenePredictionsCorrect")
	res = data.frame(as.list(res))
	return(res)
}

#### MAIN 

## read in data
genePri = fread(genePri_file, sep="\t")
colnames(genePri)[colnames(genePri)=="Disease"] = "trait" # rename Disease -> trait
traits = unique(genePri$trait)

## credible set dfs 
csPerTrait = dplyr::select(genePri, CredibleSet, trait) %>% group_by(trait) %>%  summarize(nCredibleSetsTotal = n_distinct(CredibleSet)) %>% ungroup()

uniq_cs = dplyr::select(genePri, CredibleSet, trait, TargetGene, truth) %>% dplyr::filter(truth==TRUE) %>%
	mutate(TruthGene = TargetGene) %>%
	dplyr::select(-c(TargetGene,truth)) # columns: CredibleSet, trait, TruthGene

PoPS_genes = dplyr::select(genePri, CredibleSet, trait, TargetGene, POPS.Rank) %>%
	dplyr::filter(POPS.Rank<=numPoPSGenes) %>%
	dplyr::select(-c(POPS.Rank)) %>%
	left_join(uniq_cs, by=c("CredibleSet", "trait")) # columns: CredibleSet, trait, TargetGene, TruthGene

distance_genes =  dplyr::select(genePri, CredibleSet, trait, TargetGene, PromoterDistanceToBestSNP) %>%
	arrange(trait, CredibleSet, PromoterDistanceToBestSNP) %>% group_by(CredibleSet, trait) %>%
	mutate(distanceRank = row_number()) %>%
	ungroup() %>%
	dplyr::filter(distanceRank <= numPoPSGenes) %>% # top X by distance to TSS
	dplyr::select(-c(PromoterDistanceToBestSNP,distanceRank)) %>%
	left_join(uniq_cs, by=c("CredibleSet", "trait")) # columns: CredibleSet, trait, TargetGene, TruthGene

## initialize results lists
n_rows = length(traits)
res_list_PoPS = vector("list", n_rows)
res_list_distance = vector("list", n_rows)

## iterate through traits
for (t in 1:length(traits)){
	res_list_PoPS[[t]] = get_results_vector(df=PoPS_genes, trait_this=traits[t], method="PoPS")
	res_list_distance[[t]] = get_results_vector(df=distance_genes, trait_this=traits[t], method="distanceToTSS")
}

df_res_PoPS = rbindlist(res_list_PoPS) %>% data.frame()
df_res_distance = rbindlist(res_list_distance) %>% data.frame()
df_res = rbind(df_res_PoPS, df_res_distance)
df_res = left_join(df_res, csPerTrait, by="trait") # add nCredibleSetsTotal

df_res$nCredibleSetsTotal = as.integer(df_res$nCredibleSetsTotal)
df_res$nCredibleSetsWithPredictions = as.integer(df_res$nCredibleSetsWithPredictions)
df_res$nCredibleSetGenePredictions = as.integer(df_res$nCredibleSetGenePredictions)
df_res$nCredibleSetGenePredictionsCorrect = as.integer(df_res$nCredibleSetGenePredictionsCorrect)
print(summary(df_res))

## calculate precision, recall, variance, etc.
df_res = mutate(df_res,
	recall = nCredibleSetGenePredictionsCorrect/nCredibleSetsTotal,
	precision = nCredibleSetGenePredictionsCorrect/nCredibleSetGenePredictions)

df_res = add_recall_linking_statistics(df_res, p_threshold, "baseline")
df_res = add_precision_linking_statistics(df_res, p_threshold, "baseline")

fwrite(df_res, out_file, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
