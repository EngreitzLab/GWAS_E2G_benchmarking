suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(stringr)})

main <- function() {
  option_list <- list(
    make_option(c("--input"), type="character", default=NA, help="input file"),
    make_option(c("--genes"), type="character", default=NA, help="file of genes (HGNC symbol) to filter input to; cols = chr,start,end,gene"),
    make_option(c("--biosample"), type="character", default = NA, help="biosample for this prediction file"),
    make_option(c("--invert"), type="character", default = "False", help="invert score?"),
    make_option(c("--score_col"), type="numeric", default = 6, help="col # with predictor score"))
  
  opt = parse_args(OptionParser(option_list=option_list))
  infile=opt$input;
  gene.file=opt$genes;
  biosample=opt$biosample; 
  invert=opt$invert;

  df = read.table(infile, sep="\t", header=FALSE, fill=TRUE)
  colnames(df) = c("chr", "start", "end", "TargetGene", "score")

  # read ABC data
  genes = read.table(gene.file, sep="\t", header=FALSE, fill=TRUE)
  if (ncol(genes)>1){
	colnames(genes) = c("chr", "start", "end", "hgnc.ID", "score", "strand")
  } else {
    colnames(genes)[1] = 'hgnc.ID'
  }

  # filter input file to given gene universe
  df = dplyr::filter(df, TargetGene %in% genes$hgnc.ID)
  
  # invert score if necessary
  if (invert %in% c("True", "TRUE")){
    df$score = -df$score
  }

  # set biosample name if necessary
  df$biosample = biosample
  
  # order columns
 df = df[,c('chr', 'start', 'end', 'biosample', 'TargetGene', 'score')]

  # write table
  write.table(df, file="", sep="\t", quote=F, row.names=F, col.names=F)
  }
  
  
main()
