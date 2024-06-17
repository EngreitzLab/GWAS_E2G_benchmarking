suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(colorspace)
  library(data.table)
})

main <- function() {
  ## get files from snakemake
  user.inputs = (snakemake@params$user_inputs)
  print(user.inputs)
  methods = (snakemake@params$names) %>% strsplit(" ") %>% unlist
  outFile = (snakemake@output$colorPalette)
  methods_config = fread(file=snakemake@params$methods_config, header=TRUE, sep="\t")
  
  cp = data.frame(method=methods, hex="")
  cp = left_join(cp, methods_config, by=c("method"))
  
  count = 0;
  for (i in 1:nrow(cp)){
    if (user.inputs[i]!="None"){
      cp$hex[i]=user.inputs[i]
    } else {
      count=count+1
    }
  }
  
  # generate colorspace palette and fill in 
  cols = qualitative_hcl(n=count, palette="Set 2"); 
  count=1;
  for (i in 1:nrow(cp)){
    if (cp$hex[i]==""){
      cp$hex[i] = cols[count]
      count=count+1
    }
  }
  
  # add standards
  #standards = data.frame(method=c("TSS within 100 kb", "Closest TSS", "Closest gene body", "Proximity"),
  #                       hex = c("#7A7A7A", "#A8A8A8", "#CDCDCD", "#E2E2E2"))
  #cp = rbind(cp, standards)
  
  #cpListPlotting = split(f=cp$pred_name_long, x=cp$hex)

  cp = dplyr::select(cp, method, pred_name_long, hex)
  
  # write table
  write.table(cp, outFile, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
  
}

main()
