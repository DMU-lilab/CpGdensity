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

CorPlot <- function(density, score){
  dt.cor <- data.table(density = density, score = score)
  dt.cor <- dt.cor[,.(score.median = median(score)), by = density]
  
  cor.value <- cor(dt.cor$density, dt.cor$score.median)
  
  ggplot(dt.cor, aes(x = density, y = score.median))  + geom_point() + 
    scale_x_continuous(breaks = seq(0, 250, by = 50), expand = c(0, 0), limits = c(0, 270)) +
    scale_y_continuous(breaks = seq(0, 0.8, by = 0.2), expand = c(0, 0), limits = c(0, 0.8)) +
     xlab("") + ylab("") +
    annotate("text", x = 200, y = 0.7, label = paste("cor = ", round(cor.value, 2), sep = ""), alpha = 1, fontface = "bold", size = 10) +
    publication.theme
}

load("~/Rendata/bonemarrow/mtbr_cg/chr19.Rdata")
genome <- Mmusculus
chr.name <- "chr19"
ref.fa <- genome[[chr.name]]
ref.length <- length(ref.fa)
kWinSize <- 1250

cg.mtbr <- read.csv("~/cgDensity/mm9RandomCG/random.chr19.1.mtbr")

density <- GetDensity(cg.mtbr, kWinSize, ref.length)
score <- GetScore(cg.mtbr, kWinSize, ref.length)

publication.theme <- theme_bw(15) + theme(axis.title.y=element_text(vjust=1.7), axis.title.x=element_text(vjust=-0.1), text= element_text(size = 24, face = "bold"), axis.line = element_line(colour = "black", size = 0.1), panel.border = element_rect(colour = "black", size = 1, fill = NA), panel.background = element_rect(fill = "white", size = 5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

CorPlot(density, score)
