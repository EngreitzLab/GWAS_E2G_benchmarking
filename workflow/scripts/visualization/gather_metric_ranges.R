suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
    	library(tibble)
		library(data.table)
		library(stringr)})

get_min_values <- function(df, metric_col, metric_name){
	min_val = min(df[[metric_col]])
	res =c(df$method[1], metric_name, "min", min_val)
	names(res) = c("method", "metric", "min_max", "value")
	return(as_tibble(as.list(res)))
}

get_max_values <- function(df, metric_col, metric_name){
	max_val = max(df[[metric_col]])
	res =c(df$method[1], metric_name, "max", max_val)
	names(res) = c("method", "metric", "min_max", "value")
	return(as_tibble(as.list(res)))
}

get_variant_overlap_numbers <- function(fpath, p_threshold){
	df = fread(fpath, sep="\t")
	df = dplyr::filter(df, is.finite(enrichment), is.finite(recall),!is.na(enrichment), !is.na(recall),
		nVariantsTotal>25, nCommonVariantsOverlappingEnhancers/nCommonVariantsTotal>0.001)

	res_list = vector("list", 2) # one df per metric

	for (i in 1:2){
		if (i==1){
			metric_col = "enrichment"
			metric_name = "enrichment_overlap"
			df_metric = dplyr::filter(df, p_adjust_enr<p_threshold)
		} else {
			metric_col = "recall"
			metric_name = "recall_overlap"
			df_metric = df
		}
		res_list_metric = vector("list", 2) # per min/max

		res_list_metric[[1]] = get_min_values(df_metric, metric_col, metric_name)
		res_list_metric[[2]] = get_max_values(df_metric, metric_col, metric_name)
		
		res_list[[i]] = as_tibble(rbindlist(res_list_metric))
	}

	return(as_tibble(rbindlist(res_list))) # (min, max) x (enr, recall)
}

get_gene_linking_numbers <- function(fpath, intersectPoPS){
	df = fread(fpath, sep="\t")
	df = dplyr::filter(df, is.finite(precision), is.finite(recall), !is.na(precision), !is.na(recall), nCredibleSetsTotal>0)

	metric_cols = c("precision", "recall")
	PoPS = c(TRUE, FALSE)
	res_list = vector("list", length(metric_cols)*length(PoPS)) # one df per metric, plus or minus PoPS, 
	ind = 1
	for (i in 1:length(metric_cols)){
		metric_col = metric_cols[i]
		for (k in 1:length(PoPS)){
			intersectPoPS = PoPS[k]
			if (intersectPoPS){
				metric_name = paste0(metric_col, "_linking_intersectPoPS")
				df_metric = dplyr::filter(df, intersectPoPS)
			} else {
				metric_name = paste0(metric_col, "_", "linking")
				df_metric = dplyr::filter(df, !intersectPoPS)
			}
			res_list_metric = vector("list", 2) # per min/max
			res_list_metric[[1]] = get_min_values(df_metric, metric_col, metric_name)
			res_list_metric[[2]] =  get_max_values(df_metric, metric_col, metric_name)
 
			res_list[[ind]] = as_tibble(rbindlist(res_list_metric))
			ind = ind + 1
		}
	}

	return(as_tibble(rbindlist(res_list))) # (min, max) x (prec, recall) x (with and without PoPS)
}


main <- function() {
	## process data 
	varOvl_files = snakemake@input$variantOverlap %>% strsplit(" ") %>% unlist()
	geneLinking_files = snakemake@input$geneLinking %>% strsplit(" ") %>% unlist()
	p_threshold = snakemake@params$p_threshold %>% as.numeric()
	out_file = snakemake@output$ranges

	## create df with columns: method, metric, min_max, value
	n_methods = length(varOvl_files)
	res_list = vector("list", n_methods*2)
	for (i in 1:n_methods){
		res_list[[i*2-1]]= get_variant_overlap_numbers(fpath = varOvl_files[i], p_threshold)
		res_list[[i*2]] = get_gene_linking_numbers(fpath = geneLinking_files[i])
	}

	df_res = as_tibble(rbindlist(res_list))
	fwrite(df_res, out_file, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}

main()
