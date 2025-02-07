suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  source(snakemake@params$helpful_math)})

## files from snakemake
TSS_file = snakemake@params$TSS
varInt_files = snakemake@input$variantsInt %>% strsplit(" ") %>% unlist()
varIntGroup_files = snakemake@input$variantsIntGroups %>% strsplit(" ") %>% unlist()
genePri_file = snakemake@params$genePrioritizationTable
biosamples = snakemake@params$biosamples %>% strsplit(" ") %>% unlist()
groups = snakemake@params$biosampleGroups %>% strsplit(" ") %>% unlist()
p_threshold = snakemake@params$thresholdPval %>% as.numeric()
numPoPSGenes = snakemake@params$numPoPSGenes %>% as.numeric()
numPredGenes = snakemake@params$numPredGenes %>% as.numeric()
boolean = (snakemake@params$boolean)
method_this = snakemake@wildcards$method
out_file = snakemake@output$res

#### HELPER FUNCTIONS
# return a df row with results & correct column names
get_results_vector  <- function(df, biosample, group, trait, intersectPoPS){
	nCredibleSetsOverlappingEnhancers = dplyr::select(df, CredibleSet, TruthGene) %>% distinct() %>% nrow() %>% as.integer()
	nCredibleSetsOverlappingEnhancersAnyGene = dplyr::select(df, CredibleSet, TargetGene) %>% distinct() %>% nrow() %>% as.integer()
	nCredibleSetsOverlappingEnhancersCorrectGene = dplyr::select(df, CredibleSet, TargetGene, TruthGene) %>%
		dplyr::filter(TargetGene==TruthGene) %>% distinct() %>% nrow() %>% as.integer()

	res = c(trait, biosample, group, intersectPoPS, nCredibleSetsOverlappingEnhancers, nCredibleSetsOverlappingEnhancersAnyGene, nCredibleSetsOverlappingEnhancersCorrectGene)
	names(res) = c("trait", "biosample", "group", "intersectPoPS", "nCredibleSetsOverlappingEnhancers", "nCredibleSetsOverlappingEnhancersAnyGene", "nCredibleSetsOverlappingEnhancersCorrectGene")
	res = as_tibble(as.list(res))
	return(res)
}

#### MAIN 

## read in data
TSS <- fread(TSS_file, col.names = c("chr", "start", "end", "name", "score", "strand", "Ensembl_ID", "gene_type"))
genePri = fread(genePri_file, sep="\t") %>%
	filter(TargetGene %in% TSS$name)
genePri = fread(genePri_file, sep="\t")
colnames(genePri)[colnames(genePri)=="Disease"] = "trait" # rename Disease -> trait
traits = unique(genePri$trait)

## create dataframe with the following columns: biosample[/group], trait, intersectPoPS(bool), nCredibleSets(perTrait), nCredibleSetsOverlappingEnhancers, nCredibleSetsOverlappingEnhancersAnyGene, nCredibleSetsOverlappingEnhancersCorrectGene, group(bool)
# get number of credible sets per trait and a data table of unique credible sets
csPerTrait = dplyr::select(genePri, CredibleSet, trait) %>% group_by(trait) %>%  summarize(nCredibleSetsTotal = n_distinct(CredibleSet))
uniq_cs = dplyr::select(genePri, CredibleSet, trait, TargetGene, truth) %>%
	dplyr::filter(truth==TRUE) %>%
	mutate(TruthGene = TargetGene)%>%
	dplyr::select(-c(TargetGene,truth)) # columns: CredibleSet, trait, TruthGene

PoPS_genes = dplyr::select(genePri, CredibleSet, trait, TargetGene, POPS.Score) %>%
	group_by(CredibleSet, trait) %>% 
	arrange(desc(POPS.Score)) %>% 
	mutate(POPS.Rank = row_number()) %>% 
	ungroup() %>%
	dplyr::filter(POPS.Rank<=numPoPSGenes) %>%
	mutate(PoPSGene = TargetGene) %>%
	dplyr::select(-c(TargetGene,POPS.Rank)) # columns: CredibleSet, trait, PoPSGene

# get number of credible sets overl enhancers
# get nVariantsOverlappingEnhancers
varInt_files_all = c(varInt_files, varIntGroup_files)
biosamples_all = c(biosamples, groups)

# initialize results list
n_rows = length(traits) * length(biosamples_all) * 2
res_list = vector("list", n_rows)
ind = 1
for (b in 1:length(biosamples_all)){
	biosample_this = biosamples_all[b]
	group_this = ifelse(b>length(biosamples), TRUE, FALSE)
	varInt_biosample = fread(varInt_files_all[b], sep="\t", header=TRUE)  #(variant 1-6) chr,start,end,rsid,pip,CredibleSet,trait; (prediction 7-12) chr,start,end,biosample,TargetGene,predScore
	varInt_biosample = dplyr::filter(varInt_biosample, trait %in% traits) # remove trait groups

	for (t in 1:length(traits)){
		trait_this = traits[t]
		varInt_trait = dplyr::filter(varInt_biosample, trait == trait_this)
		if (nrow(varInt_trait)>0){
			cs_this <- dplyr::filter(uniq_cs, trait == trait_this)

			# without PoPS: filter to top two-scoring genes linked by enhancers per credible set
			pred_only = dplyr::filter(varInt_trait, CredibleSet %in% cs_this$CredibleSet) %>%
				group_by(CredibleSet, TargetGene) %>%
				summarize(maxPredScore = max(predScore)) %>% # columns = credible set / trait / score
				group_by(CredibleSet) %>%
				arrange(-maxPredScore) %>%
				#mutate(predRank = row_number()) %>%
				mutate(predRank = rank(-maxPredScore, ties.method = "min")) %>%
				ungroup() %>%
				dplyr::filter(predRank <= numPredGenes) %>%
				left_join(cs_this, by=c("CredibleSet")) # add TruthGene(s)
			print(pred_only)
			res_list[[ind]] = get_results_vector(df=pred_only, biosample=biosample_this, group=group_this, trait=trait_this, intersectPoPS=FALSE)
			ind = ind + 1

			# with PoPS: filter to top two-scoring genes linked by enhancer AND by PoPS per credible set
			pred_PoPS = left_join(pred_only, PoPS_genes, by=c("CredibleSet", "trait")) %>%  # add PoPSGene(s)
				dplyr::filter(TargetGene==PoPSGene) # only keep enhancers linked to the PoPS genes
			res_list[[ind]] = get_results_vector(df=pred_PoPS, biosample=biosample_this, group=group_this, trait=trait_this, intersectPoPS=TRUE)
			ind = ind + 1
		}
	}
}
df_res = as_tibble(rbindlist(res_list))
df_res = left_join(df_res, csPerTrait, by="trait")
df_res[is.na(df_res)] = 0

df_res$nCredibleSetsTotal = as.integer(df_res$nCredibleSetsTotal)
df_res$nCredibleSetsOverlappingEnhancers = as.integer(df_res$nCredibleSetsOverlappingEnhancers)
df_res$nCredibleSetsOverlappingEnhancersAnyGene = as.integer(df_res$nCredibleSetsOverlappingEnhancersAnyGene)
df_res$nCredibleSetsOverlappingEnhancersCorrectGene = as.integer(df_res$nCredibleSetsOverlappingEnhancersCorrectGene)

## calculate precision, recall, variance, etc.
df_res = mutate(df_res,
	recall = nCredibleSetsOverlappingEnhancersCorrectGene/nCredibleSetsTotal,
	precision = nCredibleSetsOverlappingEnhancersCorrectGene/nCredibleSetsOverlappingEnhancersAnyGene)
df_res[is.na(df_res)] = 0
df_res = add_recall_linking_statistics(df_res, p_threshold, "enhancers")
df_res = add_precision_linking_statistics(df_res, p_threshold, "enhancers")
df_res[is.na(df_res)] = 0

df_res$method = method_this

fwrite(df_res, out_file, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
