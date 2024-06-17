suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

predFiles = snakemake@input$predFiles %>% as.character() %>% strsplit(" ") %>% unlist()
geneScoreFile = snakemake@params$geneScore
outFile = snakemake@output$predictionsAggregated

# read in gene scores
geneScore = fread(geneScoreFile, header=FALSE, sep="\t")
colnames(geneScore) =c("TargetGene", "geneScore")

# iterate through predictions
for (i in 1:length(predFiles)){
	predThis = fread(predFiles[i], sep="\t", header=FALSE)
	colnames(predThis) = c("chr", "start", "end", "biosample", "TargetGene", "predScore")

	predThis = left_join(predThis, geneScore, by="TargetGene")
	predThis$geneScore = replace_na(predThis$geneScore, 0)

	if (i==1){df =predThis} else {
		df = rbind(df, predThis)
	}
}

df = dplyr::select(df, chr, start, end, geneScore)

fwrite(df, outFile, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
  


