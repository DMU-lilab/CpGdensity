library("BSgenome.Mmusculus.UCSC.mm9")
library("data.table")
library("methyutils")
library("reshape2")
library("ggplot2")


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

RescaleData <- function(density, score) {
  density.scale <- (density - min(density)) / (max(density) - min(density))
  score.scale <- (score - min(score)) / (max(score) - min(score))
  
  return(data.frame(density = density.scale, score = score.scale))
}

NormalizeData <- function(density, score) {
  # step 1, scale data to [0, 100]
  
  density.max <- density[order(density)][as.integer(length(density) * 0.995)]
  density.scale <- density * (100 / density.max)
  density.scale[density.scale > 100] <- 100 
  
  score.scale <- score * 100
  return(list(density = density.scale, score = score.scale))
}

PlotData <- function(data){
  p.data <- melt(data, id="pos") 
  y.min <- p.data$value
  y.max <- p.data$value
  y.min[p.data$variable == "density"] <- 0
  y.min[p.data$variable == "score"] <- (1 - y.min[p.data$variable == "score"])
  y.max[p.data$variable == "score"] <- 1
  p.data$y.min <- y.min
  p.data$y.max <- y.max
  ggplot(p.data, aes(x = pos)) + geom_ribbon(aes(ymin=y.min,ymax=y.max, fill = variable), alpha=0.9) + scale_fill_manual(values = c("#66b2ff", "#ff9999"))  + publication.theme + xlab("") + ylab("") + theme(legend.position="none") 
}


load("~/Rendata/bonemarrow/mtbr_cg/chr19.Rdata")

genome <- Mmusculus
ref.fa <- genome[["chr19"]]
density <- GetDensity(cg.mtbr, 1250, length(ref.fa))
score <- GetScore(cg.mtbr, 1250, length(ref.fa))

rescale.data <- RescaleData(density, score)
normalize.data <- NormalizeData(density, score)

row.plotdata <- data.frame(pos = 1:length(ref.fa), density = density, score = score)
rescale.data <- data.frame(pos = 1:length(ref.fa), density = rescale.data$density, score = rescale.data$score)
normalize.plotdata <- data.frame(pos = 1:length(ref.fa), density = normalize.data$density, score = normalize.data$score)

#################
plotdata$icom <- abs(plotdata$density + plotdata$score - 100)
plotdata$maker <- plotdata$icom < 15

region <- regionAsBed(plotdata$maker, cf.length = 1000, tolerance.length = 100,chrom = "chr19")
region$width <- region$end - region$start + 1

plot.start <- region$start[15]
plot.end <- region$end[15]
#################
plot.start <- 6044218
plot.end <- 6092778
#chr19:5,078,805-5,119,910
##chr19:5,064,361-5,123,553
#chr19:5,687,462-5,773,434
##chr19:6,044,218-6,092,778
plot.data <- data.frame(pos = c(1:(plot.end - plot.start + 1)),
                        density = normalize.plotdata$density[plot.start : plot.end], 
                        score = normalize.plotdata$score[plot.start : plot.end])
rescale.data <- RescaleData(density[plot.start:plot.end], score[plot.start:plot.end])
plot.data <- data.frame(pos = c(1:(plot.end - plot.start + 1)),
                        density = rescale.data$density, 
                        score = rescale.data$score)

PlotData(plot.data)

publication.theme <- theme_bw(15) + theme(axis.title.y=element_text(vjust=1.7), axis.title.x=element_text(vjust=-0.1), text= element_text(size = 24, face = "bold"), axis.line = element_line(colour = "black", size = 0.1), panel.border = element_rect(colour = "black", size = 1, fill = NA), panel.background = element_rect(fill = "white", size = 5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

row.plot.data <- data.frame(pos = c(1:(plot.end - plot.start + 1)),
                            density = density[plot.start : plot.end], 
                            score = score[plot.start : plot.end])
g.density <- ggplot(row.plot.data,aes(x = pos, y = density)) + geom_area(fill = "#66b2ff")  + xlab("") + ylab("") + theme(axis.title.x = element_blank())+ publication.theme
g.score <- ggplot(row.plot.data,aes(x = pos, y = score)) + geom_area( fill = "#ff9999") +  scale_y_continuous(trans="reverse")+ xlab("") + ylab("") +  publication.theme
grid.arrange(g.score, g.density)
