library(optparse)

option_list = list(
		     make_option(c("-f", "--infile"), type="character", default=NULL, 
				               help="dataset file name", metavar="character"),
		     make_option(c("-w", "--width"), type="numeric", default=1000, 
				               help="window width around TSS (bp)", metavar="numeric"),
		     make_option(c("-o", "--outputDir"), type="character", default=NULL, 
				               help="output directory", metavar="character")
		       )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$infile)){
	  print_help(opt_parser)
  stop("Please specify an input file.", call.=FALSE)
}


f <- read.table(opt$infile, h = F, sep = "\t")
sample_name <- strsplit(tail(strsplit(opt$infile, split = "/")[[1]], 1), split = "\\.")[[1]][1]

pdf(file = sprintf("%s/%s.pdf", opt$outputDir, sample_name))
# insert here the right command line to create your plot
##position / enrichissement calculé dans fichier normalized
plot(x=f$V1, y=f$V2, main=sample_name)
dev.off()
