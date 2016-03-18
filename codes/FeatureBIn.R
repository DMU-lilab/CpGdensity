##
library("BSgenome.Mmusculus.UCSC.mm9")
library("data.table")
library("methyutils")
library("ggplot2")
library("grid")

GetDensity <- function(cg.mtbr, kWinSize,ref.length) {
  colnames(cg.mtbr) <- c("chr", "posi", "rC_n", "rC_p", "rT_n", "rT_p")
  posi <- cg.mtbr$posi
  rt <- logical(ref.length)
  rt[posi] <- TRUE
  win <- list(L = kWinSize, R = kWinSize)
  return(swsCalc(rt, win))
}

GetScore <- function(cg.mtbr, kWinSize, ref.length) {
  ##mtbr score sliding windows
  colnames(cg.mtbr) <- c("chr", "posi", "rC_n", "rC_p", "rT_n", "rT_p")
  
  cg.mtbr$rC <- cg.mtbr$rC_p + cg.mtbr$rC_n
  cg.mtbr$rT <- cg.mtbr$rT_p + cg.mtbr$rT_n
  
  rC <- integer(ref.length)
  rC[cg.mtbr$posi] <- cg.mtbr$rC
  rT <- integer(ref.length)
  rT[cg.mtbr$posi] <- cg.mtbr$rT
  win <- list(L = kWinSize, R = kWinSize)
  rCs <- swsCalc(rC, win)
  rTs <- swsCalc(rT, win)
  score <- rCs/(rCs + rTs)
  score[is.na(score[])] <- 0
  
  return(score)
}

CalcBinId <- function(len, bin.size, strand = "+"){
  bin.length <- floor(len / bin.size)
  bin.remind <- len - bin.length * bin.size
  
  reverse <- (strand == "-")
  if(reverse) {
    bin.ids <- c(rep(bin.size, bin.remind),rep(bin.size:1, each = bin.length))
  } else {
    bin.ids <- c(rep(1:bin.size, each = bin.length), rep(bin.size, bin.remind))
  }
  
  return(bin.ids)
}

FeatureBin <- function(df.feature, genome, bin.size = 50, by.cg = TRUE){
  dt.feature <- data.table(df.feature)
  chr.names <- paste("chr", c(1:19), sep = "")
  bins <- NULL
  for (chr.name in chr.names){
    dt.chr.feature <- dt.feature[chrom == chr.name, ]
    chr.cg.pos <- data.frame(chrom = chr.name, pos = start(matchPattern("CG", genome[[chr.name]])))
    if(nrow(dt.chr.feature) == 0){
      next
    }
    
    if(by.cg){
      dt.chr.feature.expand <- dt.chr.feature[, .(chrom = chr.name, pos = c(start:end)), by = c("id", "strand")]
      chr.cgs <- merge(dt.chr.feature.expand, chr.cg.pos, by = c("chrom", "pos"))
      dt.chr.cgs <- data.table(chr.cgs)
      chr.bins <- dt.chr.cgs[, .(chrom = chr.name, pos = pos, bin.id = CalcBinId(length(pos), bin.size, strand)), by = id]
    } else {
      dt.chr.feature.expand <- dt.chr.feature[, .(chrom = chr.name, pos = c(start:end), bin.id = CalcBinId(end - start + 1, bin.size, strand)), by = id]
      chr.bins <- merge(dt.chr.feature.expand, chr.cg.pos, by = c("chrom", "pos"))
    }
    
    bins <- rbind(bins, chr.bins) 
  }
  
  return(bins)
}
CalcGroupId <- function(dt.feature.values, n.top = 3, n.group = 10){
  dt.density.values <- dt.feature.values[,.(chrom, density, id, bin.id)]
  dt.density.mean.values <- dt.density.values[,.(density = mean(density)), by = c("id", "chrom", "bin.id")]
  dt.density.mean.order <- dt.density.mean.values[, .(density = density[order(density,decreasing = TRUE)], bin.id = bin.id[order(density, decreasing = TRUE)]), by = c("id", "chrom")]	
  dt.density.mean.top <- dt.density.mean.order[, .(feature.mean = mean(head(density, n.top))), by = c("id", "chrom")]
  dt.density.mean <- dt.density.mean.top[order(feature.mean, decreasing = TRUE),]
  group.step <- floor(nrow(dt.density.mean) / n.group)
  group.remind <- nrow(dt.density.mean) - group.step * n.group
  dt.density.mean$group.id <- c(rep(1:n.group, each = group.step), rep(n.group, group.remind))
  return(dt.density.mean)
}


FeatureValues <- function(feature.bins, genome, win.size = 1250, mtbr.path){
  if(!file.exists(mtbr.path)){
    stop("mtbr file path \"", mtbr.path ,"\" does not exist.")
  }
  
  chr.names <- paste("chr", c(1:19), sep = "")
  feature.values <- NULL
  for (chr.name in chr.names){
    chr.feature.bins <- feature.bins[chrom == chr.name, ]
    if(nrow(chr.feature.bins) == 0){
      next
    }
    
    mtbr.filename <- paste(mtbr.path, "/", chr.name, ".Rdata", sep = "")
    if(!file.exists(mtbr.path)){
      warning("mtbr file path \"", mtbr.path ,"\" does not exist.")
      next
    }
    load(mtbr.filename)
    ref.seq <- genome[[chr.name]]
    chr.density <- GetDensity(cg.mtbr, win.size, length(ref.seq))
    chr.score <- GetScore(cg.mtbr, win.size, length(ref.seq))
    df.chr.value <- data.frame(pos = 1:length(ref.seq), density = chr.density, score = chr.score)
    chr.feature.values <- merge(df.chr.value, chr.feature.bins, by = "pos") 
    
    feature.values <- rbind(feature.values, chr.feature.values)
  }
  
  dt.feature.values <- data.table(feature.values)
  dt.feature.group <- CalcGroupId(dt.feature.values, 3, 10)
  dt.feature.result <- merge(dt.feature.values, dt.feature.group, by = c("id", "chrom"))
  dt.feature.result <- dt.feature.result[, .(density = mean(density), score = mean(score)), by = c("group.id", "bin.id")]
  return(dt.feature.result)
  
}

###############################

genome <- Mmusculus
genes <- read.csv("~/database/mm9/mm9.gene.csv",header=T,sep="\t")
genes$width <- genes$txEnd - genes$txStart
genes <- genes[genes$width > 1000,]
#genes <- genes[genes$strand == "+",]
gene.feature <- data.frame(id = genes$name, chrom = genes$chrom, start = genes$txStart, end = genes$txEnd, strand = genes$strand)

df.feature <- gene.feature

feature.bins <- FeatureBin(df.feature, genome, bin.size = 50, by.cg = TRUE)

mtbr.path <- "~/Rendata/bonemarrow/mtbr_cg/"
dt.feature.result <- FeatureValues(feature.bins, genome, win.size = 1250, mtbr.path)


dt.feature.result$group.id2 <- paste("GR",dt.feature.result$group.id,sep="")

#ggplot(dt.feature.result, aes(x = bin.id, y = score, color = group.id2)) + geom_point() + geom_line() + scale_colour_discrete(breaks = paste("GR", c(1:10), sep = ""), labels = paste("GR",c(1:10))) + theme(legend.position = c(1,1), legend.justification = c(1,0.95)) + theme(legend.background = element_blank()) + theme(legend.key = element_blank()) + publication.theme

ggplot(Line.dt.feature.result, aes(x = bin.id, y = score, color = as.factor(group.id))) + geom_point() + geom_line() + publication.theme + xlab("") + ylab("") + theme(legend.position="none") 
#write.table(dt.feature.result,"~/cgDensity/FeatureBin/dt.feature.resultchr1_19bygenelength.txt",col.names = T,row.names = F,quote = F,sep="\t")

#ggplot(dt.feature.result, aes(x = bin.id, y = density, color = group.id2)) + geom_point() + geom_line() + scale_colour_discrete(breaks = paste("GR", c(1:10), sep = ""), labels = paste("a",c(1:10))) 
 #ggplot(dt.feature.result, aes(x = bin.id, y = density, color = group.id2, group = group.id2)) + geom_point() + geom_line() 
#g <- g + theme(axis.text.x = element_text(angle = 90))
#print(g)

#library(grid)
#newtheme <- theme_bw(15) + theme(panel.border = element_rect(colour = "black", size = 0.5), panel.background = element_rect(fill = "white", size = 5), panel.grid.minor = element_blank(), axis.ticks.length = unit(0.2,"cm"))

publication.theme <- theme_bw(15) + theme(axis.title.y=element_text(vjust=1.7), axis.title.x=element_text(vjust=-0.1), text= element_text(size = 24, face = "bold"), axis.line = element_line(colour = "black", size = 0.1), panel.border = element_rect(colour = "black", size = 1, fill = NA), panel.background = element_rect(fill = "white", size = 5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

