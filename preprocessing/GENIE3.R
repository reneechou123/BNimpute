#!/usr/bin/env Rscript
library(GENIE3)
library(optparse)

option_list = list(
  make_option(c("-e", "--exp"), type="character", default=NULL,
              help="expression file name", metavar="character"),
  make_option(c("-p", "--prior"), type="character", default=NULL,
              help="prior file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

exprMatrix <- as.matrix(read.table(opt$exp, header=FALSE, row.names=1))
prior <- read.table(opt$prior, header=TRUE, row.names=1)

set.seed(123)
weightMatrix <- t(GENIE3(exprMatrix, regulators=colnames(prior)))
write.table(weightMatrix, opt$out, sep='\t', row.names=TRUE, col.names=NA, quote=FALSE)

