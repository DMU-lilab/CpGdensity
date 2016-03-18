#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("tools", quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))

# Get command line options & arguments

arguments <- parse_args(OptionParser(usage = "%prog [options] mtbr.Rdata"), positional_arguments = 1)
kMtbrFile <- arguments$args

if(!file.exists(kMtbrFile)){
  stop("mtbr file \"", kMtbrFile ,"\" does not exist.")
}

filename.base <- file_path_sans_ext(basename(kMtbrFile))
output.filename <- paste("./", filename.base, ".csv", sep = "")

load(kMtbrFile)
cg.mtbr$posi <- as.integer(cg.mtbr$posi)
cg.mtbr$rC_n <- as.integer(cg.mtbr$rC_n)
cg.mtbr$rC_p <- as.integer(cg.mtbr$rC_p)
cg.mtbr$rT_n <- as.integer(cg.mtbr$rT_n)
cg.mtbr$rT_p <- as.integer(cg.mtbr$rT_p)

write.csv(cg.mtbr, output.filename, row.names = FALSE, quote = FALSE)
