library(optparse)

option_list = list(
		     make_option(c("-f", "--infile"), type="character", default=NULL, 
				               help="dataset file name", metavar="character"),
		     make_option(c("-o", "--outputDir"), type="character", default=NULL, 
				               help="output directory", metavar="character")
		       )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$infile)){
	  print_help(opt_parser)
  stop("Please specify an input file.", call.=FALSE)
}




f <- read.table(opt$infile, h=F)
head(f)
sum(f$V1)
sample_name <- strsplit(tail(strsplit(opt$infile, split = "/")[[1]], 1), split = "\\.")[[1]][1]

#counts <- rep(f$V2, f$V1)

pdf(file = sprintf("%s/%s.pdf", opt$outputDir, sample_name))
#plot(x = f$V2, y = f$V1, main = sample_name, xlab = "Insert size", ylab = "n", xlim = c(0,1000))
#plot(density(counts, bw = 0.5))
plot(density(f$V1, bw = 0.05), xlim=c(0,600))
#on change echelle pour voir quelque chose

dev.off()
