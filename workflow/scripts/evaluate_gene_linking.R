suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(GenomicRanges)
  library(R.utils)})

# get info from snakemake
pred_file = snakemake@input$predictionsMerged
PoPS_file = snakemake@input$benchmarkingData
traits_use = snakemake@params$traitsToUse %>% as.character() %>% strsplit(" ") %>% unlist()
bgVariant_files = snakemake@input$bgVariantFiles %>% as.character() %>% strsplit(" ") %>% unlist()
bgVariant_chroms = snakemake@params$bgVariantChroms %>% as.character() %>% str_replace("chr", "") %>% strsplit(" ") %>% unlist() # numbers only
nBootstrap = snakemake@params$nBootstrap %>% as.numeric()
UKBB_traits_file = snakemake@params$UKBB_SuSiE_traits
UKBB_cs_dir = snakemake@params$UKBB_SuSiE_csDir
numPoPSGenes = snakemake@wildcards$numPoPSGenes %>% as.numeric()
numPredGenes = snakemake@wildcards$numPredGenes %>% as.numeric()
usePoPS = snakemake@wildcards$usePoPS %>% as.character()
tissue = snakemake@wildcards$tissue
method = snakemake@wildcards$method
outFile = snakemake@output$outputTable

## READ IN DATA
eg_preds = fread(pred_file, sep="\t", header=TRUE)

# fine-mapped variant set
traits = fread(UKBB_traits_file, header=FALSE) %>% setNames("trait")
traits = traits$trait
if (!("ALL" %in% traits_use)){ # filter if necessary
	traits = intersect(traits, traits_use)
} 
n_traits = length(traits)
traitVariantsList = vector(mode="list", n_traits)
for (i in 1:n_traits){
	fname = file.path(UKBB_cs_dir, traits[i], "variant.list.txt")
	temp = fread(fname, sep="\t", header=TRUE)
	traitVariantsList[[i]] = temp
}
names(traitVariantsList) = traits

# background variants
n_chr = length(bgVariant_chroms)
bgVarList = vector(mode="list", n_chr)
for (i in 1:n_chr){
	temp = fread(bgVariant_files[i], sep="\t", header=FALSE)
	colnames(temp) = c("CHR", "SNP", "CM", "BP", "REF", "ALT")
	bgVarList[[i]] = temp
}
names(bgVarList) = paste0("chr", bgVariant_chroms)

# PoPS data
pops_silver = fread(PoPS_file, sep="\t", header=TRUE)
if (!("ALL" %in% traits_use)){
	pops_silver = dplyr::filter(pops_silver, Disease %in% traits_use)
}

## OVERLAP PREDICTIONS WITH CAUSAL-GENE-LINKED REGIONS
message("before gr1")
gr1 = GRanges(seqnames = eg_preds$chr, ranges=IRanges(start=eg_preds$start, end=eg_preds$end))

chrs = as.character(sapply(pops_silver$CredibleSet, function(x) return(strsplit(x, ":")[[1]][1])))
cred_start =  as.numeric(sapply(pops_silver$CredibleSet, function(x) return(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])))
cred_end =  as.numeric(sapply(pops_silver$CredibleSet, function(x) return(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2])))

pops_silver$chrs = chrs
pops_silver$cred_start = cred_start
pops_silver$cred_end = cred_end

#ucreds = unique(pops_silver$CredibleSet)
ucreds = dplyr::select(pops_silver, Disease, CredibleSet) %>% distinct()
message("Number of unique CS-trait pairs: ", nrow(ucreds))

true_positives = c()
positives = c()
relevant = c()
#creds_vec = c() # not used

# iterate through credible sets
for(numu in 1:nrow(ucreds)){
#for (numu in 1:10){ # for testing
	this_cs = ucreds$CredibleSet[numu]; this_trait = ucreds$Disease[numu]
	temp = dplyr::filter(pops_silver, CredibleSet==this_cs, Disease==this_trait)
	causal_gene = temp$TargetGene[which(temp$truth == TRUE)] # can be multiple for different diseases
    #trait_names = unique(temp$Disease) # not used anywhere
	#print(trait_names)
	#print(ucreds[numu])

    if(usePoPS=="withPoPS"){
		pops_genes = temp$TargetGene[order(temp$POPS.Rank, decreasing = F)[1:numPoPSGenes]] 
	}

	message("before gr_cred")
    chr_cred = as.character(strsplit(this_cs, ":")[[1]][1])
    # rsids = c()
    # pips = c()

	# # iterate through traits - actually iterate through diesase_name and get PoPS genes within loop?
    # for(num_trait in 1:length(trait_names)){
	# 	this_vars = traitVariantsList[[trait_names[num_trait]]]
	# 	idx = which(this_vars$CredibleSet==ucreds[numu])
	# 	message("Trait: ", trait_names[num_trait], ", CS: ", ucreds[numu], ", Length idx: ", length(idx))
	# 	rsids = c(rsids, this_vars$rsid[idx])
	# 	pips = c(pips, this_vars$pip[idx])
    # }

	this_vars = traitVariantsList[[this_trait]]
	idx = which(this_vars$CredibleSet == this_cs)
	message("Trait: ", this_trait, ", CS: ", this_cs, ", Length idx: ", length(idx))
	rsids = this_vars$rsid[idx]
	pips = this_vars$pip[idx]

	message("# all rsids: ", length(rsids))
    rsids = rsids[which(pips > 0.10)]
	message("# rsids with PIP>0.1: ", length(rsids))
	message("# intersecting with bg variants: ", length(intersect(rsids,bgVarList[[chr_cred]]$SNP)))
    #coding_promoter_snps = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Promoter_Coding_variants.txt", header=F)[,1]
    #rsids = setdiff(rsids, coding_promoter_snps)

    if(length(intersect(rsids,bgVarList[[chr_cred]]$SNP)) == 0){
		message("no intersection")
    	next
    }
    #starts = bgVarList[[chr_cred]]$BP[match(intersect(rsids,bgVarList[[chr_cred]]$SNP), bgVarList[[chr_cred]]$SNP)]-1
    #ends = bgVarList[[chr_cred]]$BP[match(intersect(rsids,bgVarList[[chr_cred]]$SNP), bgVarList[[chr_cred]]$SNP)]+1

	bgVar_this = dplyr::filter(bgVarList[[chr_cred]], SNP %in% rsids)
	print(bgVar_this)

    gr_cred = GRanges(seqnames = chr_cred, ranges = IRanges(start = bgVar_this$BP-1, end = bgVar_this$BP+1))
    cc=findOverlaps(gr_cred, gr1, type = "any", select = "all")

    if(length(cc) > 0){
      any = c(any, 1)
      eg_preds3 = eg_preds[unique(subjectHits(cc)), c("TargetGene", "predScore")]
      gene_scores = tapply(eg_preds3$predScore, eg_preds3$TargetGene, max)
      gene_scores2 = rep(0, length(temp$TargetGene))
      names(gene_scores2) = temp$TargetGene
      gene_scores2[match(intersect(names(gene_scores), temp$TargetGene), temp$TargetGene)] = gene_scores[match(intersect(names(gene_scores), temp$TargetGene), names(gene_scores))]

      if(max(gene_scores2) == 0){
        true_positives = c(true_positives, 0)
        positives = c(positives, 0)
      }else{
        genes_selected = names(gene_scores2)[order(gene_scores2, decreasing = T)[1:pmin(numPredGenes, length(gene_scores2))]]
        if(usePoPS=="withPoPS"){
          genes_selected = intersect(genes_selected, pops_genes)
        }
        true_positives = c(true_positives, length(intersect(genes_selected, causal_gene)))
        positives = c(positives, length(genes_selected))
      }
    }else{
      any = c(any, 0)
      true_positives = c(true_positives, 0)
      positives = c(positives, 0)
    }
    relevant = c(relevant, length(causal_gene))
    cat("Processing credible set:", numu, "out of ", nrow(ucreds), "\n")
  }

## PRECISION AND RECALL
precision = sum(true_positives[which(positives !=  0)])/length(positives[which(positives !=  0)]) # length positives
recall = sum(true_positives[which(relevant !=  0)])/length(relevant[which(relevant !=  0)])

precision2 = sum(true_positives[which(positives !=  0)])/sum(positives[which(positives !=  0)]) # sum positives, use this one
recall2 = sum(true_positives[which(relevant !=  0)])/sum(relevant[which(relevant !=  0)])

precision_boot = c()
precision2_boot = c()
recall_boot = c()
recall2_boot = c()

for(nboot in 1:nBootstrap){
   idx = sample(1:length(positives), length(positives), replace = T)
   true_positives_boot = true_positives[idx]
   positives_boot = positives[idx]
   relevant_boot =- relevant[idx]

   precision_boot = c(precision_boot, sum(true_positives_boot[which(positives_boot !=  0)])/length(positives_boot[which(positives_boot !=  0)]))
   recall_boot = c(recall_boot, sum(true_positives_boot[which(relevant_boot !=  0)])/length(relevant_boot[which(relevant_boot !=  0)]))

   precision2_boot = c(precision2_boot, sum(true_positives_boot[which(positives_boot !=  0)])/sum(positives_boot[which(positives_boot !=  0)]))
   recall2_boot = c(recall2_boot, sum(true_positives_boot[which(relevant_boot !=  0)])/sum(relevant_boot[which(relevant_boot !=  0)]))
 }

sd_precision = sd(precision_boot)
sd_recall = sd(recall_boot)
sd_precision2 = sd(precision2_boot)
sd_recall2 = sd(recall2_boot)

df_final = data.frame(precision=precision, sd_precision=sd_precision, recall=recall, sd_recall=sd_recall,
									precision2=precision2, sd_precision2=sd_precision2, recall2=recall2, sd_recall2=sd_recall2,
									tissue=tissue, use_pops=usePoPS, num_pops_genes=numPoPSGenes, num_pred_genes=numPredGenes,
									method=method, num_relevant=sum(relevant))

fwrite(df_final, outFile, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
