suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(R.utils)})

annotFiles = snakemake@input$bgVariantAnnotations %>% as.character() %>% strsplit(" ") %>% unlist()
traitFiles = snakemake@params$finemappedVarFiles %>% as.character() %>% strsplit(" ") %>% unlist()
traitNames = snakemake@params$finemappedVarTraits%>% as.character() %>% strsplit(" ") %>% unlist()
outFile = snakemake@output$results
n_boot = snakemake@params$nBootstrap %>% as.numeric()
method = snakemake@wildcards$method

## read in files
n_chr = length(annotFiles)
annotList = vector(mode = "list", n_chr) # pred file per chrom
for (i in 1:n_chr){
	temp = fread(annotFiles[i], sep="\t", header=TRUE)
	annotList[[i]] = temp
}

n_traits = length(traitFiles)
trait_var_count = data.frame(trait=traitNames, num_variants=0)
traitList = vector(mode="list", n_traits) #variant file per trait
for (i in 1:n_traits){
	temp = fread(traitFiles[i], header=FALSE)[[1]] # list of rsids
	traitList[[i]] = temp
	trait_var_count$num_variants[i] = length(temp)
}

## initialize some variables
sum_farh_chr = rep(0, n_traits)
sum_finemap = rep(0, n_traits)
sum_binannot = 0
sum_all = 0
sum_aggr_annot = matrix(0, n_boot, n_traits)
per_chr_1 = matrix(0, n_chr, n_traits) # chrom by trait matrix
per_chr_2 = array(0, n_chr)
per_chr_3 = matrix(0, n_chr, n_traits)

## STEP ONE: compile data across variant files + chromosomes
for (chr in 1:n_chr){	# iterate through chromosomes
	df_annot = annotList[[chr]] # "df1" 
	binannot = df_annot$ANNOT # vector of annotations for each prediction
	all_trait_annotations = vector(mode="list", n_traits)

	for (tr in 1:n_traits) { 	# iterate through traits
		df_trait = unlist(traitList[[tr]]) # "finemap1" 
		print(head(df_annot))
		print(head(df_trait))
		
		df_annot$this_annot <- ifelse(df_annot$SNP %in% df_trait, 1, 0)
		this_annot = df_annot$this_annot
		#this_annot = rep(0, length(binannot)) # 0s for each element in this chr, "finemap1_annot"
		#print(head(intersect(df_annot$SNP, df_trait)))
		#print(head(match(intersect(df_annot$SNP, df_trait), df_annot$SNP)))
		#this_annot[match(intersect(df_annot$SNP, df_trait), df_annot$SNP)] = 1 
		print(sum(this_annot))
		all_trait_annotations[[tr]] = this_annot # save for bootstrapping
		per_chr_1[chr, tr] = sum(binannot[this_annot==1])
		per_chr_3[chr, tr] = sum(this_annot)
	}

	# for all traits in a chromosomes
	sum_farh_chr = sum_farh_chr + per_chr_1[chr,]
	per_chr_2[chr] = sum(binannot) # s																			
	sum_binannot = sum_binannot + per_chr_2[chr]
	sum_all = sum_all + nrow(df_annot) # number of bg snps overlapping pred
	sum_finemap = sum_finemap + per_chr_3[chr,]

	# bootstrapping 
	aggr_annot = c()
	for (n in 1:n_boot){ # boostrap
		binannot_samp = sample(binannot, length(binannot), replace = F)
		this_row = rep(0, n_traits)
		for (tr in 1:n_traits){ # SO MANY LOOPS
			this_annot_again = all_trait_annotations[[tr]]
			this_row[tr] = sum(binannot_samp[this_annot_again==1])
		}
		aggr_annot = rbind(aggr_annot, this_row)
	}
	sum_aggr = sum_aggr_annot + aggr_annot	
	message("We are at chromosome ", chr, "/", n_chr)
}

## NEXT STEP: enrichment, precision, recall
enr_bootmat=matrix(0, n_boot, n_traits)
prec_bootmat=matrix(0, n_boot, n_traits)
recall_bootmat=matrix(0, n_boot, n_traits)

for (nb in 1:n_boot){
	idx = sample(1:n_chr, n_chr, replace=TRUE) # permutation of chromosomes 
	print(idx)
	print(per_chr_1[idx, ])

	enr_bootmat[nb, ] = (colSums(per_chr_1[idx, ])/colSums(per_chr_3[idx, ]))/(sum(per_chr_2[idx])/sum_all)
  	prec_bootmat[nb, ] = (colSums(per_chr_1[idx, ])/sum(per_chr_2[idx]))
	recall_bootmat[nb, ] = (colSums(per_chr_1[idx, ])/colSums(per_chr_3[idx, ]))
}

enr = (sum_farh_chr/sum_binannot)/(sum_finemap/sum_all)
precision = (sum_farh_chr/sum_binannot)
recall = (sum_farh_chr/sum_finemap)

senr = apply(enr_bootmat, 2, sd)
sprec = apply(prec_bootmat, 2, sd)
srecall = apply(recall_bootmat, 2, sd)


boot_enr = matrix(0, n_boot, n_traits)
boot_prec = matrix(0, n_boot, n_traits)
boot_recall = matrix(0, n_boot, n_traits)

for (mm  in 1:n_boot){
  boot_enr[mm, ] = (sum_aggr_annot[mm, ]/sum_binannot)/(sum_finemap/sum_all)
  boot_prec[mm, ] = (sum_aggr_annot[mm, ]/sum_binannot)
  boot_recall[mm, ] = (sum_aggr_annot[mm, ]/sum_finemap)
}

penr = c()
pprec = c()
precall = c()

for(mm in 1:n_traits){
  penr = c(penr, pnorm(enr[mm], mean(boot_enr[,mm]), sd(boot_enr[,mm]), lower.tail = F))
  pprec = c(pprec, pnorm(precision[mm], mean(boot_prec[,mm]), sd(boot_prec[,mm]), lower.tail = F))
  precall = c(precall, pnorm(recall[mm], mean(boot_recall[,mm]), sd(boot_recall[,mm]), lower.tail = F))
}

df_final = data.frame(cbind(traitNames, enr, senr, penr, precision, sprec, pprec, recall, srecall, precall))
colnames(df_final) =  c("trait", "enrichment", "sd_enrichment", "pval_enrichment", "precision", "sd_precision", "pval_precision", "recall", "sd_recall", "pval_recall")
df_final$method = method
df_final = left_join(df_final, trait_var_count, by="trait") 

fwrite(df_final, outFile, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
  

#################################################################

