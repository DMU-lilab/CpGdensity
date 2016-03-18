library("ggplot2")
library("tools")
library("BSgenome.Mmusculus.UCSC.mm9")
library("IRanges")
library("GenomicRanges")

source("~/work/cgdensity/cgdensity.lib.R")

CountCmpl <- function(ref.filename, mtbr.filename, sw.size = 1250, cmpl.threshhold = 0.3) {
  ref.seq <- readDNAStringSet(ref.filename)[[1]]
  cg.pos <- start(matchPattern("CG", ref.seq))
  mtbr <- read.csv(mtbr.filename)
  ref.length <- length(ref.seq)
  
  cg.density <- GetDensity(mtbr, ref.length, sw.size)
  cg.score <- GetScore(mtbr, ref.length, sw.size)
  
  df.all.rescale <- RescaleData(cg.density, cg.score)
  df.cg.rescale <- df.all.rescale[cg.pos, ]
  
  v.cmpl <- df.cg.rescale$density + df.cg.rescale$score - 1
  count.cmpl <- sum(v.cmpl > -cmpl.threshhold & v.cmpl <= 0)
  count.uncmpl <- sum(v.cmpl > 0 | v.cmpl < -cmpl.threshhold)
  
  return(data.frame(count.cmpl = count.cmpl, count.uncmpl = count.uncmpl))
}

GetAllCmpl <- function(file.path, output.filename) {
  if(!file.exists(file.path)) {
    stop(file.path, " does not exist!")
  }
  
  file.names <- list.files(file.path, "*.fa")
  file.bases <- file_path_sans_ext(file.names)
  df.cmpl <- NULL
  for(file.base in file.bases) {
    ref.filename <- paste(file.path, "/", file.base, ".fa", sep = "")
    mtbr.filename <- paste(file.path, "/", file.base, ".mtbr", sep = "")
    
    if(!file.exists(ref.filename)) {
      warning(ref.filename, " does not exist!")
      next
    }
    
    if(!file.exists(mtbr.filename)) {
      warning(mtbr.filename, " does not exist!")
      next
    }
    
    message("processing ", file.base, " ", date())
    cmpl <- CountCmpl(ref.filename, mtbr.filename, 1250, 0.3)
    cmpl$sample <- file.base
    df.cmpl <- rbind(df.cmpl, cmpl)
    
    write.table(cmpl, file = output.filename, append = TRUE, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
  }
  
  return(df.cmpl)
}

random.data.path <- "~/work/cgdensity/data/bonemarrow/random.data/random.chr19/"
df.cmpl <- GetAllCmpl(random.data.path, "./cmpl.count.csv")