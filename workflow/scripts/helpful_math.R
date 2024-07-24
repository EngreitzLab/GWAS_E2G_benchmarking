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

get_n_equals <- function(df, type){
	if (type=="variant_overlap"){
		col_name = "nVariantsTotal"
	} else if (type=="gene_linking") {
		col_name = "nCredibleSetsTotal"
	}

	n_min = min(df[[col_name]])
	n_max = max(df[[col_name]])
	if (n_min==n_max){
			n_label = paste0("N = ", n_min)
	} else {
		n_label = paste0("N = ", n_min, "-", n_max)
	}
}
