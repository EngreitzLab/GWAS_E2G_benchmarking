suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

# CI + pval (hypergeometric) for enrichment aka risk ratio > 1
# CI ref: https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_confidence_intervals/bs704_confidence_intervals8.html
add_enrichment_statistics <- function(df, p_threshold){
	z = qnorm(p_threshold/2, lower.tail=FALSE) # e.g. 1.96

	df = mutate(df, 
		SE_log_enr  = sqrt(((nVariantsTotal-nVariantsOverlappingEnhancers)/nVariantsOverlappingEnhancers)/nVariantsTotal) + ((nCommonVariantsTotal-nCommonVariantsOverlappingEnhancers)/nCommonVariantsOverlappingEnhancers)/nCommonVariantsTotal,
		CI_enr_low = exp(log(enrichment)-z*SE_log_enr),
		CI_enr_high = exp(log(enrichment)+z*SE_log_enr),
		p_enr = phyper(nVariantsOverlappingEnhancers, nVariantsTotal, nCommonVariantsTotal, nVariantsOverlappingEnhancers + nCommonVariantsOverlappingEnhancers, log.p=FALSE, lower.tail=FALSE),
		p_adjust_enr = p.adjust(p_enr, method="bonferroni"))

		return(df)
}

# CI for recall (Wilson score interval): https://users.stat.ufl.edu/~aa/articles/agresti_coull_1998.pdf
add_recall_overlap_statistics <- function(df, p_threshold){
	z = qnorm(p_threshold/2, lower.tail=FALSE) # e.g. 1.96

	df = mutate(df,
		recall_adjust = (nVariantsOverlappingEnhancers+2)/(nVariantsTotal+4),
		SE_recall = sqrt(recall_adjust * (1-recall_adjust)/(nVariantsTotal+4)),
		CI_recall_low = recall - z*SE_recall,
		CI_recall_high = recall + z*SE_recall)

	return(df)
}

summarize_grouped_enrichment_recall <- function(df, p_threshold){
	z = qnorm(p_threshold/2, lower.tail=FALSE) # e.g. 1.96

	df = summarize(df,
			nVariantsTotal = sum(nVariantsTotal),
			nVariantsOverlappingEnhancers = sum(nVariantsOverlappingEnhancers),
			nCommonVariantsTotal = sum(nCommonVariantsTotal),
			nCommonVariantsOverlappingEnhancers = sum(nCommonVariantsOverlappingEnhancers)) %>%
		mutate(recall = nVariantsOverlappingEnhancers/nVariantsTotal,
			enrichment = recall/(nCommonVariantsOverlappingEnhancers/nCommonVariantsTotal))

	df = add_enrichment_statistics(df, p_threshold)
	df = add_recall_overlap_statistics(df, p_threshold)
	
	return(df)
}

# based on Wilson score interval as above
add_recall_linking_statistics <- function(df, p_threshold, type){
	z = qnorm(p_threshold/2, lower.tail=FALSE) # e.g. 1.96

	if (type=="enhancers"){
		TP_col = "nCredibleSetsOverlappingEnhancersCorrectGene"
	} else if (type=="baseline") {
		TP_col = "nCredibleSetGenePredictionsCorrect"
	}

	df = df %>% mutate(
		recall_adjust = (!!sym(TP_col)+2)/(nCredibleSetsTotal+4),
		SE_recall = sqrt(recall_adjust * (1-recall_adjust)/(nCredibleSetsTotal+4)),
		CI_recall_low = recall - z*SE_recall,
		CI_recall_high = recall + z*SE_recall)

		return(df)
}

# based on Wilson score interval as above
add_precision_linking_statistics <- function(df, p_threshold, type){
	z = qnorm(p_threshold/2, lower.tail=FALSE) # e.g. 1.96

	if (type=="enhancers"){
		TP_col = "nCredibleSetsOverlappingEnhancersCorrectGene"
		predP_col = "nCredibleSetsOverlappingEnhancersAnyGene"
	} else if (type=="baseline"){
		TP_col = "nCredibleSetGenePredictionsCorrect"
		predP_col = "nCredibleSetGenePredictions"
	}

	df = mutate(df,
		precision_adjust = (!!sym(TP_col)+2)/(!!sym(predP_col)+4),
		SE_precision = sqrt(precision_adjust * (1-precision_adjust)/(!!sym(predP_col)+4)),
		CI_precision_low = precision - z*SE_precision,
		CI_precision_high = precision + z*SE_precision)

		return(df)
}

summarize_grouped_precision_recall<- function(df, p_threshold, type){
	z = qnorm(p_threshold/2, lower.tail=FALSE) # e.g. 1.96

	if (type=="enhancers"){
		TP_col = "nCredibleSetsOverlappingEnhancersCorrectGene"
		predP_col = "nCredibleSetsOverlappingEnhancersAnyGene"
	} else if (type=="baseline"){
		TP_col = "nCredibleSetGenePredictionsCorrect"
		predP_col = "nCredibleSetGenePredictions"
	}

	df = summarize(df,
			nCredibleSetsTotal = sum(nCredibleSetsTotal),
			temp_TP = sum(!!sym(TP_col)),
			temp_predP = sum(!!sym(predP_col))) %>%
		mutate(recall = temp_TP/nCredibleSetsTotal,
			precision = temp_TP/temp_predP)
	df[[TP_col]] = df[["temp_TP"]]
	df[[predP_col]] = df[["temp_predP"]]
	df = dplyr::select(df, -c(temp_TP, temp_predP))

	df = add_recall_linking_statistics(df, p_threshold, type)
	df = add_precision_linking_statistics(df, p_threshold, type)
	
	return(df)
}

get_min_max_from_column <- function(df, col_name) {
	n_min = min(df[[col_name]])
	n_max = max(df[[col_name]])
	if (n_min==n_max){
			n_label = paste0(n_min)
	} else {
		n_label = paste0( n_min, "-", n_max)
	}

	return(n_label)
}

get_n_equals_variant_overlap <- function(df, var, traits, trait_groups) {
	# var = all variants file
	
	traits_to_use <- c()
	for (i in seq_along(traits)) {
		this_trait <- traits[i]
		if (this_trait %in% colnames(trait_groups)){
			traits_to_use <- c(traits_to_use, trait_groups[[this_trait]])
		}
		else {
			traits_to_use <- c(traits_to_use, this_trait)
		}
	}
	traits_to_use <- traits_to_use[!is.na(traits_to_use)] %>% unique()

	#chr,start,end,rsid,CredibleSet,trait 
	n_var <- dplyr::filter(var, trait %in% traits_to_use) %>% dplyr::select(chr, start, end) %>% distinct() %>% nrow()
	n_var_text <- paste0(n_var, " variants from ", length(traits_to_use), " traits")

	print(df)
	n_pairs <- dplyr::select(df, biosample, trait, nVariantsTotal) %>% distinct()
	n_pairs_text <-  paste0("(", sum(n_pairs$nVariantsTotal), " variant-biosample pairs)")

	label <- paste0(n_var_text, "\n", n_pairs_text)
	return(label)
}

get_n_equals_gene_linking <- function(df){
	# number of credible sets & traits
	n_cs <- df %>% dplyr::select(trait, nCredibleSetsTotal) %>% distinct()
	n_cs_text <- paste0(sum(n_cs$nCredibleSetsTotal), " credible sets from ", length(unique(n_cs$trait)), " traits tested")

	n_pairs <- df %>% dplyr::select(biosample, trait, nCredibleSetsTotal) %>% distinct()
	n_pairs_text <- paste0("(", sum(n_pairs$nCredibleSetsTotal), " credible set-biosample pairs)")

	label <- paste0(n_cs_text, "\n", n_pairs_text)
	return(label)
}