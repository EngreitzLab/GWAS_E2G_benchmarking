suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)})

predFiles = snakemake@input$predFiles %>% as.character() %>% strsplit(" ") %>% unlist()
outFile = snakemake@output$predictionsAggregated

# iterate through predictions
predList = vector(mode="list", length(predFiles))
for (i in 1:length(predFiles)){
	predThis = fread(predFiles[i], sep="\t", header=FALSE)
	colnames(predThis) = c("chr", "start", "end", "biosample", "TargetGene", "predScore")
	predList[[i]] = predThis
}

df =rbindlist(predList)

fwrite(df, outFile, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
  


