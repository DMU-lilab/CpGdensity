library("ggplot2")
library("zoo")

##sliding windows
SlidingWindow<-function(x, win=list(L=75, R=75)){
  xLen <- length(x)
  winLen<-win$L + win$R + 1
  sws <- rollsum(x,winLen)
  sws_head<-tail(cumsum(x[1:(winLen-1)]),win$L)
  sws_tail<-rev(tail(cumsum(rev(x[(xLen - winLen + 2):xLen])),win$R))
  sws <- append(sws_head, sws)
  sws <- append(sws, sws_tail)
  return(sws)
}

GetDensity <- function(chr.mtbr, chr.length, sw.size) {
  cgs <- chr.mtbr$pos
  v.cgs <- integer(chr.length)
  v.cgs[cgs] <- 1
  
  return (SlidingWindow(v.cgs, list(L = sw.size, R = sw.size)))
}

GetScore <- function(chr.mtbr, chr.length, sw.size) {
  rc <- with(chr.mtbr, rC_n + rC_p)
  sum <- with(chr.mtbr, rC_n + rC_p + rT_n + rT_p)
  
  v.rc <- integer(chr.length)
  v.sum <- integer(chr.length)
  
  v.rc[chr.mtbr$pos] <- rc
  v.sum[chr.mtbr$pos] <- sum
  
  v.rc <- SlidingWindow(v.rc, list(L = sw.size, R = sw.size))
  v.sum <- SlidingWindow(v.sum, list(L = sw.size, R = sw.size))
  v.score <- v.rc / v.sum
  v.score[is.na(v.score)] <- 0
  
  return (v.score)
}

NormlizeData <- function(density, score, win.size) {
  density.scale <- (density * 1000) / (win.size)
  score.scale <- score
  
  return(data.frame(density = density.scale, score = score.scale))
}

RescaleData <- function(density, score) {
  # Find outliers using Hampel procedure
  
  #density.outlier <- hampel(density, 20, 3)
  #density.scale <- density / (max(density.outlier$y) - min(density.outlier$y))
  #density.scale[density.outlier$ind] <- 1
  
  #density.outlier <- quantile(density, 0.9995)
  #density.scale <- density / density.outlier
  #density.scale[density.scale > 1] <- 1
  
  #density.outlier <- median(density) + 15 * mad(density)
  #density.scale <- density / density.outlier
  #density.scale[density.scale > 1] <- 1
  
  #score.outlier <- median(score) + 8 * mad(score)
  #score.outlier <- quantile(score, 0.9995)
  #score.scale <- score / score.outlier
  #score.scale[score.scale > 1] <- 1
  
  density.scale <- (density - min(density)) / (max(density) - min(density))
  score.scale <- (score - min(score)) / (max(score) - min(score))
  
  return(data.frame(density = density.scale, score = score.scale))
}

SliceCount <- function(density, score, range = list(start = -1.0, end = 1.0), delta = 0.02) {
  v.sum <- density + score
  v.dist <- seq(range$start, range$end, by = delta)
  v.cdf <- unlist(lapply(v.dist, function(dist) if(dist <= 0) { sum((v.sum > 1 + dist) & (v.sum <= 1))} else {sum((v.sum < 1 + dist) & (v.sum >= 1))})) 
  v.pdf <- unlist(lapply(v.dist, function(dist, delta = 0.02) sum((v.sum > 1 + dist) & (v.sum < 1 + dist + delta)), delta = delta))
  
  return(data.frame(dist = v.dist, pdf = v.pdf, cdf = v.cdf))
}
