suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
})

### INPUTS
varInt_files = snakemake@input$varInt %>% as.character() %>% strsplit(" ") %>% unlist() # (variant 1-6) chr,start,end,rsid,CredibleSet,trait; (prediction 7-12) chr,start,end,biosample,TargetGene,predScore
biosamples = snakemake@params$biosamples %>% as.character() %>% strsplit(" ") %>% unlist()
traits = snakemake@params$traits %>% as.character() %>% strsplit(" ") %>% unlist()
provided_threshold = snakemake@params$threshold %>% as.numeric()
n_steps_target = snakemake@params$nSteps %>% as.numeric()
comparisons = snakemake@params$comparisonTable
boolean = snakemake@params$boolean
out_file = snakemake@output$thresholdSpan

# for 0/1 predictors
if (boolean %in% c("TRUE", "True", TRUE)){
  print("boolean")
  n_steps_target=2
} else {
  print("non-boolean")
} 

map = data.frame(biosample=biosamples, trait=traits)
map = mutate(map, key = paste0(biosample, ".", trait))

res_list = vector("list",  length(varInt_files))
for (i in 1:length(varInt_files)){
	df = fread(varInt_files[i], sep="\t")
	res_list[[i]] = mutate(df, key = paste0(biosample, ".", trait)) %>%
		dplyr::filter(key %in% map$key)
}
res = data.frame(rbindlist(res_list)) %>%
	dplyr::select(varChr, varStart, varEnd, predScore) %>%
	dplyr::filter(!is.na(predScore), is.finite(predScore)) %>%
	group_by(-c(predScore)) %>%
	summarize(predScore = max(predScore)) %>%
	distinct()

# evenly-spaced scores (n_steps_target)
scores = as.numeric(res$predScore)
n_distinct_scores = unique(scores) %>% length()
even_spaced_scores = data.frame(threshold = seq(min(scores), max(scores), length.out=n_steps_target))

# quantile spaced scores (n_steps_target)
prob_vector = seq(0, 1, length.out=n_steps_target)
quantile_scores = unname(quantile(scores, na.rm=T, probs=prob_vector)) 
quantile_scores = c(quantile_scores, provided_threshold)

thresholds = data.frame(threshold=quantile_scores) %>% distinct()

n_steps = n_steps_target
while ((nrow(thresholds)<n_steps_target) & (n_steps<n_distinct_scores)){ # less than target # thresholds, but not max-ed out
	n_steps = n_steps + 1
	prob_vector = seq(0, 1, length.out = n_steps)
	quantile_scores = unname(quantile(scores, na.rm=T, probs=prob_vector))
	quantile_scores = c(quantile_scores, provided_threshold)

	thresholds = data.frame(threshold=quantile_scores) %>% distinct()
}

# combine
thresholds_all = rbind(thresholds, even_spaced_scores) %>%
	distinct() %>%
	arrange(threshold)

print(thresholds_all)
write.table(thresholds_all, out_file, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
