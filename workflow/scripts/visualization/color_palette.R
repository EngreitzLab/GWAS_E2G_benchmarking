suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(colorspace)
  library(data.table)
})

main <- function() {
	## get files from snakemake
	methods = (snakemake@params$names) %>% strsplit(" ") %>% unlist
	outFile = (snakemake@output$colorPalette)
	methods_config = fread(file=snakemake@params$methods_config, header=TRUE, sep="\t") %>%
		dplyr::select(method, pred_name_long, color)

	cp = data.frame(method=methods)
	cp = left_join(cp, methods_config, by=c("method")) # to get pred_name_long
	cp$hex = "#000000"

	for (i in 1:nrow(cp)){
		cp$hex[i] = ifelse(is.na(cp$color[i]), NA, cp$color[i])
	}
	num_na = sum(is.na(cp$color))
  
  # generate colorspace palette and fill in 
  cols = qualitative_hcl(n=num_na, palette="Set 2"); 
  count=1;
  for (i in 1:nrow(cp)){
    if (is.na(cp$hex[i])){
      cp$hex[i] = cols[count]
      count=count+1
    }
  }

  # add baselines
  baselines = data.frame(method=c("PoPS", "distanceToTSS"),
                         				hex = c("#1c2a43", "#c5cad7"), # dark and light grey
										pred_name_long = c("Polygenic priority score (PoPS)", "Distance to TSS"))   

  cp = dplyr::select(cp, method, pred_name_long, hex)
  cp = rbind(cp, baselines)
  
  # write table
  write.table(cp, outFile, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
  
}

main()
