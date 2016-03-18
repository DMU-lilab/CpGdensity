## Fig1 plot CMD
## ten species CG and AT Hilbert plot and histgram
library(data.table)
library(HilbertCurve)
library(circlize)
library(ComplexHeatmap)
library(methyutils)
library(ggplot2)

## Function

patternHistogramPlot <- function(chr.seq,level=9,pattern="CG",species="hg38",max.count=150,adjust=10){
  require(HilbertCurve)
  chr.seq <- chr.seq
  hilbert.level = level
  target.pattern <- pattern
  pt <- matchPattern(target.pattern,chr.seq)
  pt.pos <- start(pt)
  
  hbc = HilbertCurve(1, length(chr.seq), level = hilbert.level, mode = "pixel")
  chr.brks <- hbc@BINS@start
  pt.nu <- table(cut(pt.pos,c(chr.brks,length(chr.seq))))
  pt.bin.nu <- as.integer(pt.nu)/hbc@BINS@width*1000

  hg <- data.frame(id=c(1:length(pt.bin.nu)),val=pt.bin.nu)
  theme_set(theme_bw(12))
  p <- ggplot(hg,aes(val))+geom_histogram(fill="white", colour="black",binwidth=1)+xlim(0,max.count) + xlab(paste(species,pattern,"count"))
  return(p)
}

patternDensityPlot <- function(chr.seq,level=9,pattern="CG",species="hg38",max.count=150,adjust=10){
  require(HilbertCurve)
  chr.seq <- chr.seq
  hilbert.level = level
  target.pattern <- pattern
  pt <- matchPattern(target.pattern,chr.seq)
  pt.pos <- start(pt)
  
  hbc = HilbertCurve(1, length(chr.seq), level = hilbert.level, mode = "pixel")
  chr.brks <- hbc@BINS@start
  pt.nu <- table(cut(pt.pos,c(chr.brks,length(chr.seq))))
  pt.bin.nu <- as.integer(pt.nu)/hbc@BINS@width*1000
  
  hg <- data.frame(id=c(1:length(pt.bin.nu)),val=pt.bin.nu)
  theme_set(theme_bw(12))
  p <- ggplot(hg)+geom_density(aes(x=val),fill="RoyalBlue",adjust=adjust)+xlim(0,max.count) + xlab(paste(species,pattern,"count"))
  return(p)
}
##
wigPlot <- function(chr.seq,pattern="CG",wig.file="finalPlot/win2500.wig",win=1250,chrom="chr1"){
  l <- numeric(length(chr.seq))
  l[start(matchPattern(pattern,chr.seq))] <- 1
  cgw <- swsCalc(l ,win=list(L=win,R=win))
  head <- paste("variableStep chrom=",chrom," span=1\n", sep = "")
  win2500 <- data.frame(pos=c(1:length(chr.seq)),count=cgw)
  win2500 <- win2500[which(win2500$count != 0),]
  cat(head,file = wig.file,append=F)
  write.table(win2500,file = wig.file, row.names=F, col.names=F, quote=F, append=T)
}
##
binCount <- function(chr.seq,level=9,pattern="CG"){
  pt <- matchPattern(pattern,chr.seq)
  pt.pos <- start(pt)
  hbc = HilbertCurve(1, length(chr.seq), level = level, mode = "pixel")
  chr.brks <- hbc@BINS@start
  pt.nu <- table(cut(pt.pos,c(chr.brks,length(chr.seq))))
  pt.bin.nu <- as.integer(pt.nu)/hbc@BINS@width*1000
  return(pt.bin.nu)
}

############################

##hg38
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38

chr="chr1"
hilbert.level <- 9
chr.seq <- hg38[[chr]]
## adjust=1
hg38.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"CG count density", "adjust=1" ))
hg38.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"AT count density", "adjust=1" ))
hg38.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"GC count density", "adjust=1" ))
hg38.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"TA count density", "adjust=1" ))
hg38.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"AC count density", "adjust=1" ))
## adjust=5
hg38.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"CG count density", "adjust=5" ))
hg38.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"AT count density", "adjust=5" ))
hg38.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"GC count density", "adjust=5" ))
hg38.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"TA count density", "adjust=5" ))
hg38.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"AC count density", "adjust=5" ))
## adjust=10
hg38.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"CG count density", "adjust=10"))
hg38.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"AT count density", "adjust=10"))
hg38.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"GC count density", "adjust=10"))
hg38.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"TA count density", "adjust=10"))
hg38.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/hg38-",chr,"-pattern-count.density.pdf",sep = ""))
  
  hg38.density.CG.ad1
  hg38.density.AT.ad1
  hg38.density.GC.ad1
  hg38.density.TA.ad1
  hg38.density.AC.ad1
  
  hg38.density.CG.ad5
  hg38.density.AT.ad5
  hg38.density.GC.ad5
  hg38.density.TA.ad5
  hg38.density.AC.ad5
  
  hg38.density.CG.ad10
  hg38.density.AT.ad10
  hg38.density.GC.ad10
  hg38.density.TA.ad10
  hg38.density.AC.ad10

dev.off()

chr="chr2"
hilbert.level <- 9
chr.seq <- hg38[[chr]]
## adjust=1
hg38.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"CG count density", "adjust=1" ))
hg38.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"AT count density", "adjust=1" ))
hg38.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"GC count density", "adjust=1" ))
hg38.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"TA count density", "adjust=1" ))
hg38.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"AC count density", "adjust=1" ))
## adjust=5
hg38.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"CG count density", "adjust=5" ))
hg38.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"AT count density", "adjust=5" ))
hg38.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"GC count density", "adjust=5" ))
hg38.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"TA count density", "adjust=5" ))
hg38.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"AC count density", "adjust=5" ))
## adjust=10
hg38.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"CG count density", "adjust=10"))
hg38.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"AT count density", "adjust=10"))
hg38.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"GC count density", "adjust=10"))
hg38.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"TA count density", "adjust=10"))
hg38.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/hg38-",chr,"-pattern-count.density.pdf",sep = ""))

hg38.density.CG.ad1
hg38.density.AT.ad1
hg38.density.GC.ad1
hg38.density.TA.ad1
hg38.density.AC.ad1

hg38.density.CG.ad5
hg38.density.AT.ad5
hg38.density.GC.ad5
hg38.density.TA.ad5
hg38.density.AC.ad5

hg38.density.CG.ad10
hg38.density.AT.ad10
hg38.density.GC.ad10
hg38.density.TA.ad10
hg38.density.AC.ad10

dev.off()

chr="chr3"
hilbert.level <- 9
chr.seq <- hg38[[chr]]
## adjust=1
hg38.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"CG count density", "adjust=1" ))
hg38.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"AT count density", "adjust=1" ))
hg38.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"GC count density", "adjust=1" ))
hg38.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"TA count density", "adjust=1" ))
hg38.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"AC count density", "adjust=1" ))
## adjust=5
hg38.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"CG count density", "adjust=5" ))
hg38.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"AT count density", "adjust=5" ))
hg38.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"GC count density", "adjust=5" ))
hg38.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"TA count density", "adjust=5" ))
hg38.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"AC count density", "adjust=5" ))
## adjust=10
hg38.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"CG count density", "adjust=10"))
hg38.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"AT count density", "adjust=10"))
hg38.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"GC count density", "adjust=10"))
hg38.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"TA count density", "adjust=10"))
hg38.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/hg38-",chr,"-pattern-count.density.pdf",sep = ""))

hg38.density.CG.ad1
hg38.density.AT.ad1
hg38.density.GC.ad1
hg38.density.TA.ad1
hg38.density.AC.ad1

hg38.density.CG.ad5
hg38.density.AT.ad5
hg38.density.GC.ad5
hg38.density.TA.ad5
hg38.density.AC.ad5

hg38.density.CG.ad10
hg38.density.AT.ad10
hg38.density.GC.ad10
hg38.density.TA.ad10
hg38.density.AC.ad10

dev.off()

chr="chr4"
hilbert.level <- 9
chr.seq <- hg38[[chr]]
## adjust=1
hg38.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"CG count density", "adjust=1" ))
hg38.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"AT count density", "adjust=1" ))
hg38.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"GC count density", "adjust=1" ))
hg38.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"TA count density", "adjust=1" ))
hg38.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"AC count density", "adjust=1" ))
## adjust=5
hg38.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"CG count density", "adjust=5" ))
hg38.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"AT count density", "adjust=5" ))
hg38.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"GC count density", "adjust=5" ))
hg38.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"TA count density", "adjust=5" ))
hg38.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"AC count density", "adjust=5" ))
## adjust=10
hg38.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"CG count density", "adjust=10"))
hg38.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"AT count density", "adjust=10"))
hg38.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"GC count density", "adjust=10"))
hg38.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"TA count density", "adjust=10"))
hg38.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/hg38-",chr,"-pattern-count.density.pdf",sep = ""))

hg38.density.CG.ad1
hg38.density.AT.ad1
hg38.density.GC.ad1
hg38.density.TA.ad1
hg38.density.AC.ad1

hg38.density.CG.ad5
hg38.density.AT.ad5
hg38.density.GC.ad5
hg38.density.TA.ad5
hg38.density.AC.ad5

hg38.density.CG.ad10
hg38.density.AT.ad10
hg38.density.GC.ad10
hg38.density.TA.ad10
hg38.density.AC.ad10

dev.off()

chr="chr5"
hilbert.level <- 9
chr.seq <- hg38[[chr]]
## adjust=1
hg38.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"CG count density", "adjust=1" ))
hg38.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"AT count density", "adjust=1" ))
hg38.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"GC count density", "adjust=1" ))
hg38.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"TA count density", "adjust=1" ))
hg38.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"AC count density", "adjust=1" ))
## adjust=5
hg38.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"CG count density", "adjust=5" ))
hg38.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"AT count density", "adjust=5" ))
hg38.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"GC count density", "adjust=5" ))
hg38.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"TA count density", "adjust=5" ))
hg38.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"AC count density", "adjust=5" ))
## adjust=10
hg38.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"CG count density", "adjust=10"))
hg38.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"AT count density", "adjust=10"))
hg38.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"GC count density", "adjust=10"))
hg38.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"TA count density", "adjust=10"))
hg38.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/hg38-",chr,"-pattern-count.density.pdf",sep = ""))

hg38.density.CG.ad1
hg38.density.AT.ad1
hg38.density.GC.ad1
hg38.density.TA.ad1
hg38.density.AC.ad1

hg38.density.CG.ad5
hg38.density.AT.ad5
hg38.density.GC.ad5
hg38.density.TA.ad5
hg38.density.AC.ad5

hg38.density.CG.ad10
hg38.density.AT.ad10
hg38.density.GC.ad10
hg38.density.TA.ad10
hg38.density.AC.ad10

dev.off()

chr="chr19"
hilbert.level <- 8
chr.seq <- hg38[[chr]]
## adjust=1
hg38.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"CG count density", "adjust=1" ))
hg38.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"AT count density", "adjust=1" ))
hg38.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"GC count density", "adjust=1" ))
hg38.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"TA count density", "adjust=1" ))
hg38.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 1 ) + labs(title= paste("hg38",chr,"AC count density", "adjust=1" ))
## adjust=5
hg38.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"CG count density", "adjust=5" ))
hg38.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"AT count density", "adjust=5" ))
hg38.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"GC count density", "adjust=5" ))
hg38.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"TA count density", "adjust=5" ))
hg38.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 5 ) + labs(title= paste("hg38",chr,"AC count density", "adjust=5" ))
## adjust=10
hg38.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"CG count density", "adjust=10"))
hg38.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"AT count density", "adjust=10"))
hg38.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"GC count density", "adjust=10"))
hg38.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"TA count density", "adjust=10"))
hg38.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "hg38", max.count = 150,adjust = 10) + labs(title= paste("hg38",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/hg38-",chr,"-pattern-count.density.pdf",sep = ""))

hg38.density.CG.ad1
hg38.density.AT.ad1
hg38.density.GC.ad1
hg38.density.TA.ad1
hg38.density.AC.ad1

hg38.density.CG.ad5
hg38.density.AT.ad5
hg38.density.GC.ad5
hg38.density.TA.ad5
hg38.density.AC.ad5

hg38.density.CG.ad10
hg38.density.AT.ad10
hg38.density.GC.ad10
hg38.density.TA.ad10
hg38.density.AC.ad10

dev.off()


##mm10
library(BSgenome.Mmusculus.UCSC.mm10)
mm10 <- BSgenome.Mmusculus.UCSC.mm10

chr="chr1"
hilbert.level <- 9
chr.seq <- mm10[[chr]]
## adjust=1
mm10.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"CG count density", "adjust=1" ))
mm10.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"AT count density", "adjust=1" ))
mm10.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"GC count density", "adjust=1" ))
mm10.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"TA count density", "adjust=1" ))
mm10.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"AC count density", "adjust=1" ))
## adjust=5
mm10.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"CG count density", "adjust=5" ))
mm10.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"AT count density", "adjust=5" ))
mm10.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"GC count density", "adjust=5" ))
mm10.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"TA count density", "adjust=5" ))
mm10.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"AC count density", "adjust=5" ))
## adjust=10
mm10.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"CG count density", "adjust=10"))
mm10.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"AT count density", "adjust=10"))
mm10.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"GC count density", "adjust=10"))
mm10.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"TA count density", "adjust=10"))
mm10.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/mm10-",chr,"-pattern-count.density.pdf",sep = ""))

mm10.density.CG.ad1
mm10.density.AT.ad1
mm10.density.GC.ad1
mm10.density.TA.ad1
mm10.density.AC.ad1

mm10.density.CG.ad5
mm10.density.AT.ad5
mm10.density.GC.ad5
mm10.density.TA.ad5
mm10.density.AC.ad5

mm10.density.CG.ad10
mm10.density.AT.ad10
mm10.density.GC.ad10
mm10.density.TA.ad10
mm10.density.AC.ad10

dev.off()

chr="chr2"
hilbert.level <- 9
chr.seq <- mm10[[chr]]
## adjust=1
mm10.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"CG count density", "adjust=1" ))
mm10.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"AT count density", "adjust=1" ))
mm10.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"GC count density", "adjust=1" ))
mm10.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"TA count density", "adjust=1" ))
mm10.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"AC count density", "adjust=1" ))
## adjust=5
mm10.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"CG count density", "adjust=5" ))
mm10.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"AT count density", "adjust=5" ))
mm10.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"GC count density", "adjust=5" ))
mm10.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"TA count density", "adjust=5" ))
mm10.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"AC count density", "adjust=5" ))
## adjust=10
mm10.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"CG count density", "adjust=10"))
mm10.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"AT count density", "adjust=10"))
mm10.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"GC count density", "adjust=10"))
mm10.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"TA count density", "adjust=10"))
mm10.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/mm10-",chr,"-pattern-count.density.pdf",sep = ""))

mm10.density.CG.ad1
mm10.density.AT.ad1
mm10.density.GC.ad1
mm10.density.TA.ad1
mm10.density.AC.ad1

mm10.density.CG.ad5
mm10.density.AT.ad5
mm10.density.GC.ad5
mm10.density.TA.ad5
mm10.density.AC.ad5

mm10.density.CG.ad10
mm10.density.AT.ad10
mm10.density.GC.ad10
mm10.density.TA.ad10
mm10.density.AC.ad10

dev.off()

chr="chr3"
hilbert.level <- 9
chr.seq <- mm10[[chr]]
## adjust=1
mm10.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"CG count density", "adjust=1" ))
mm10.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"AT count density", "adjust=1" ))
mm10.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"GC count density", "adjust=1" ))
mm10.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"TA count density", "adjust=1" ))
mm10.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"AC count density", "adjust=1" ))
## adjust=5
mm10.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"CG count density", "adjust=5" ))
mm10.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"AT count density", "adjust=5" ))
mm10.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"GC count density", "adjust=5" ))
mm10.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"TA count density", "adjust=5" ))
mm10.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"AC count density", "adjust=5" ))
## adjust=10
mm10.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"CG count density", "adjust=10"))
mm10.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"AT count density", "adjust=10"))
mm10.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"GC count density", "adjust=10"))
mm10.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"TA count density", "adjust=10"))
mm10.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/mm10-",chr,"-pattern-count.density.pdf",sep = ""))

mm10.density.CG.ad1
mm10.density.AT.ad1
mm10.density.GC.ad1
mm10.density.TA.ad1
mm10.density.AC.ad1

mm10.density.CG.ad5
mm10.density.AT.ad5
mm10.density.GC.ad5
mm10.density.TA.ad5
mm10.density.AC.ad5

mm10.density.CG.ad10
mm10.density.AT.ad10
mm10.density.GC.ad10
mm10.density.TA.ad10
mm10.density.AC.ad10

dev.off()

chr="chr4"
hilbert.level <- 9
chr.seq <- mm10[[chr]]
## adjust=1
mm10.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"CG count density", "adjust=1" ))
mm10.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"AT count density", "adjust=1" ))
mm10.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"GC count density", "adjust=1" ))
mm10.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"TA count density", "adjust=1" ))
mm10.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"AC count density", "adjust=1" ))
## adjust=5
mm10.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"CG count density", "adjust=5" ))
mm10.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"AT count density", "adjust=5" ))
mm10.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"GC count density", "adjust=5" ))
mm10.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"TA count density", "adjust=5" ))
mm10.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"AC count density", "adjust=5" ))
## adjust=10
mm10.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"CG count density", "adjust=10"))
mm10.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"AT count density", "adjust=10"))
mm10.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"GC count density", "adjust=10"))
mm10.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"TA count density", "adjust=10"))
mm10.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/mm10-",chr,"-pattern-count.density.pdf",sep = ""))

mm10.density.CG.ad1
mm10.density.AT.ad1
mm10.density.GC.ad1
mm10.density.TA.ad1
mm10.density.AC.ad1

mm10.density.CG.ad5
mm10.density.AT.ad5
mm10.density.GC.ad5
mm10.density.TA.ad5
mm10.density.AC.ad5

mm10.density.CG.ad10
mm10.density.AT.ad10
mm10.density.GC.ad10
mm10.density.TA.ad10
mm10.density.AC.ad10

dev.off()

chr="chr5"
hilbert.level <- 9
chr.seq <- mm10[[chr]]
## adjust=1
mm10.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"CG count density", "adjust=1" ))
mm10.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"AT count density", "adjust=1" ))
mm10.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"GC count density", "adjust=1" ))
mm10.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"TA count density", "adjust=1" ))
mm10.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"AC count density", "adjust=1" ))
## adjust=5
mm10.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"CG count density", "adjust=5" ))
mm10.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"AT count density", "adjust=5" ))
mm10.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"GC count density", "adjust=5" ))
mm10.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"TA count density", "adjust=5" ))
mm10.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"AC count density", "adjust=5" ))
## adjust=10
mm10.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"CG count density", "adjust=10"))
mm10.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"AT count density", "adjust=10"))
mm10.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"GC count density", "adjust=10"))
mm10.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"TA count density", "adjust=10"))
mm10.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/mm10-",chr,"-pattern-count.density.pdf",sep = ""))

mm10.density.CG.ad1
mm10.density.AT.ad1
mm10.density.GC.ad1
mm10.density.TA.ad1
mm10.density.AC.ad1

mm10.density.CG.ad5
mm10.density.AT.ad5
mm10.density.GC.ad5
mm10.density.TA.ad5
mm10.density.AC.ad5

mm10.density.CG.ad10
mm10.density.AT.ad10
mm10.density.GC.ad10
mm10.density.TA.ad10
mm10.density.AC.ad10

dev.off()

chr="chr19"
hilbert.level <- 8
chr.seq <- mm10[[chr]]
## adjust=1
mm10.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"CG count density", "adjust=1" ))
mm10.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"AT count density", "adjust=1" ))
mm10.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"GC count density", "adjust=1" ))
mm10.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"TA count density", "adjust=1" ))
mm10.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 1 ) + labs(title= paste("mm10",chr,"AC count density", "adjust=1" ))
## adjust=5
mm10.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"CG count density", "adjust=5" ))
mm10.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"AT count density", "adjust=5" ))
mm10.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"GC count density", "adjust=5" ))
mm10.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"TA count density", "adjust=5" ))
mm10.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 5 ) + labs(title= paste("mm10",chr,"AC count density", "adjust=5" ))
## adjust=10
mm10.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"CG count density", "adjust=10"))
mm10.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"AT count density", "adjust=10"))
mm10.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"GC count density", "adjust=10"))
mm10.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"TA count density", "adjust=10"))
mm10.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "mm10", max.count = 150,adjust = 10) + labs(title= paste("mm10",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/mm10-",chr,"-pattern-count.density.pdf",sep = ""))

mm10.density.CG.ad1
mm10.density.AT.ad1
mm10.density.GC.ad1
mm10.density.TA.ad1
mm10.density.AC.ad1

mm10.density.CG.ad5
mm10.density.AT.ad5
mm10.density.GC.ad5
mm10.density.TA.ad5
mm10.density.AC.ad5

mm10.density.CG.ad10
mm10.density.AT.ad10
mm10.density.GC.ad10
mm10.density.TA.ad10
mm10.density.AC.ad10

dev.off()

## taeGut1
library(BSgenome.Tguttata.UCSC.taeGut1)
taeGut1 <- BSgenome.Tguttata.UCSC.taeGut1

chr="chr1"
hilbert.level <- 8
chr.seq <- taeGut1[[chr]]
## adjust=1
taeGut1.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=1" ))
taeGut1.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=1" ))
taeGut1.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=1" ))
taeGut1.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=1" ))
taeGut1.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=1" ))
## adjust=5
taeGut1.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=5" ))
taeGut1.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=5" ))
taeGut1.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=5" ))
taeGut1.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=5" ))
taeGut1.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=5" ))
## adjust=10
taeGut1.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=10"))
taeGut1.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=10"))
taeGut1.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=10"))
taeGut1.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=10"))
taeGut1.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/taeGut1-",chr,"-pattern-count.density.pdf",sep = ""))

  taeGut1.density.CG.ad1
  taeGut1.density.AT.ad1
  taeGut1.density.GC.ad1
  taeGut1.density.TA.ad1
  taeGut1.density.AC.ad1
  
  taeGut1.density.CG.ad5
  taeGut1.density.AT.ad5
  taeGut1.density.GC.ad5
  taeGut1.density.TA.ad5
  taeGut1.density.AC.ad5
  
  taeGut1.density.CG.ad10
  taeGut1.density.AT.ad10
  taeGut1.density.GC.ad10
  taeGut1.density.TA.ad10
  taeGut1.density.AC.ad10

dev.off()

chr="chr2"
hilbert.level <- 8
chr.seq <- taeGut1[[chr]]
## adjust=1
taeGut1.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=1" ))
taeGut1.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=1" ))
taeGut1.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=1" ))
taeGut1.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=1" ))
taeGut1.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=1" ))
## adjust=5
taeGut1.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=5" ))
taeGut1.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=5" ))
taeGut1.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=5" ))
taeGut1.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=5" ))
taeGut1.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=5" ))
## adjust=10
taeGut1.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=10"))
taeGut1.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=10"))
taeGut1.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=10"))
taeGut1.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=10"))
taeGut1.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/taeGut1-",chr,"-pattern-count.density.pdf",sep = ""))

  taeGut1.density.CG.ad1
  taeGut1.density.AT.ad1
  taeGut1.density.GC.ad1
  taeGut1.density.TA.ad1
  taeGut1.density.AC.ad1
  
  taeGut1.density.CG.ad5
  taeGut1.density.AT.ad5
  taeGut1.density.GC.ad5
  taeGut1.density.TA.ad5
  taeGut1.density.AC.ad5
  
  taeGut1.density.CG.ad10
  taeGut1.density.AT.ad10
  taeGut1.density.GC.ad10
  taeGut1.density.TA.ad10
  taeGut1.density.AC.ad10

dev.off()

chr="chr3"
hilbert.level <- 8
chr.seq <- taeGut1[[chr]]
## adjust=1
taeGut1.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=1" ))
taeGut1.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=1" ))
taeGut1.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=1" ))
taeGut1.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=1" ))
taeGut1.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=1" ))
## adjust=5
taeGut1.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=5" ))
taeGut1.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=5" ))
taeGut1.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=5" ))
taeGut1.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=5" ))
taeGut1.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=5" ))
## adjust=10
taeGut1.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=10"))
taeGut1.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=10"))
taeGut1.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=10"))
taeGut1.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=10"))
taeGut1.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/taeGut1-",chr,"-pattern-count.density.pdf",sep = ""))

  taeGut1.density.CG.ad1
  taeGut1.density.AT.ad1
  taeGut1.density.GC.ad1
  taeGut1.density.TA.ad1
  taeGut1.density.AC.ad1
  
  taeGut1.density.CG.ad5
  taeGut1.density.AT.ad5
  taeGut1.density.GC.ad5
  taeGut1.density.TA.ad5
  taeGut1.density.AC.ad5
  
  taeGut1.density.CG.ad10
  taeGut1.density.AT.ad10
  taeGut1.density.GC.ad10
  taeGut1.density.TA.ad10
  taeGut1.density.AC.ad10

dev.off()

chr="chr4"
hilbert.level <- 8
chr.seq <- taeGut1[[chr]]
## adjust=1
taeGut1.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=1" ))
taeGut1.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=1" ))
taeGut1.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=1" ))
taeGut1.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=1" ))
taeGut1.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=1" ))
## adjust=5
taeGut1.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=5" ))
taeGut1.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=5" ))
taeGut1.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=5" ))
taeGut1.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=5" ))
taeGut1.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=5" ))
## adjust=10
taeGut1.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=10"))
taeGut1.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=10"))
taeGut1.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=10"))
taeGut1.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=10"))
taeGut1.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/taeGut1-",chr,"-pattern-count.density.pdf",sep = ""))

taeGut1.density.CG.ad1
taeGut1.density.AT.ad1
taeGut1.density.GC.ad1
taeGut1.density.TA.ad1
taeGut1.density.AC.ad1

taeGut1.density.CG.ad5
taeGut1.density.AT.ad5
taeGut1.density.GC.ad5
taeGut1.density.TA.ad5
taeGut1.density.AC.ad5

taeGut1.density.CG.ad10
taeGut1.density.AT.ad10
taeGut1.density.GC.ad10
taeGut1.density.TA.ad10
taeGut1.density.AC.ad10

dev.off()

chr="chr5"
hilbert.level <- 8
chr.seq <- taeGut1[[chr]]
## adjust=1
taeGut1.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=1" ))
taeGut1.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=1" ))
taeGut1.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=1" ))
taeGut1.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=1" ))
taeGut1.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=1" ))
## adjust=5
taeGut1.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=5" ))
taeGut1.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=5" ))
taeGut1.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=5" ))
taeGut1.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=5" ))
taeGut1.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=5" ))
## adjust=10
taeGut1.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=10"))
taeGut1.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=10"))
taeGut1.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=10"))
taeGut1.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=10"))
taeGut1.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/taeGut1-",chr,"-pattern-count.density.pdf",sep = ""))

taeGut1.density.CG.ad1
taeGut1.density.AT.ad1
taeGut1.density.GC.ad1
taeGut1.density.TA.ad1
taeGut1.density.AC.ad1

taeGut1.density.CG.ad5
taeGut1.density.AT.ad5
taeGut1.density.GC.ad5
taeGut1.density.TA.ad5
taeGut1.density.AC.ad5

taeGut1.density.CG.ad10
taeGut1.density.AT.ad10
taeGut1.density.GC.ad10
taeGut1.density.TA.ad10
taeGut1.density.AC.ad10

dev.off()

chr="chr19"
hilbert.level <- 8
chr.seq <- taeGut1[[chr]]
## adjust=1
taeGut1.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=1" ))
taeGut1.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=1" ))
taeGut1.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=1" ))
taeGut1.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=1" ))
taeGut1.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 1 ) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=1" ))
## adjust=5
taeGut1.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=5" ))
taeGut1.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=5" ))
taeGut1.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=5" ))
taeGut1.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=5" ))
taeGut1.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 5 ) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=5" ))
## adjust=10
taeGut1.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"CG count density", "adjust=10"))
taeGut1.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"AT count density", "adjust=10"))
taeGut1.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"GC count density", "adjust=10"))
taeGut1.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"TA count density", "adjust=10"))
taeGut1.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "taeGut1", max.count = 150,adjust = 10) + labs(title= paste("taeGut1",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/taeGut1-",chr,"-pattern-count.density.pdf",sep = ""))

taeGut1.density.CG.ad1
taeGut1.density.AT.ad1
taeGut1.density.GC.ad1
taeGut1.density.TA.ad1
taeGut1.density.AC.ad1

taeGut1.density.CG.ad5
taeGut1.density.AT.ad5
taeGut1.density.GC.ad5
taeGut1.density.TA.ad5
taeGut1.density.AC.ad5

taeGut1.density.CG.ad10
taeGut1.density.AT.ad10
taeGut1.density.GC.ad10
taeGut1.density.TA.ad10
taeGut1.density.AC.ad10

dev.off()

## parkeri
library("Biostrings")
parkeri.path <- "~/projects/CpGden/parkeri.fasta"
parkeri <- readDNAStringSet(parkeri.path, format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=F)
hilbert.level <- 9
pkr <- DNAString()
for(tg.seq in c(1:100) ){
  # pkr <- DNAString(c(pkr,parkeri[tg.seq]))
  pkr <- DNAStringSet(paste(as.character(pkr),as.character(parkeri[tg.seq]),sep=""))
}
chr.seq <- pkr[[1]]

parkeri.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "parkeri", max.count = 150,adjust = 1 ) + labs(title= paste("parkeri top100 CG count density", "adjust=1" ))
parkeri.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "parkeri", max.count = 150,adjust = 1 ) + labs(title= paste("parkeri top100 CG count density", "adjust=1" ))
parkeri.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "parkeri", max.count = 150,adjust = 1 ) + labs(title= paste("parkeri top100 CG count density", "adjust=1" ))
parkeri.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "parkeri", max.count = 150,adjust = 1 ) + labs(title= paste("parkeri top100 CG count density", "adjust=1" ))
parkeri.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "parkeri", max.count = 150,adjust = 1 ) + labs(title= paste("parkeri top100 CG count density", "adjust=1" ))

parkeri.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "parkeri", max.count = 150,adjust = 5 ) + labs(title= paste("parkeri top100 CG count density", "adjust=5" ))
parkeri.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "parkeri", max.count = 150,adjust = 5 ) + labs(title= paste("parkeri top100 CG count density", "adjust=5" ))
parkeri.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "parkeri", max.count = 150,adjust = 5 ) + labs(title= paste("parkeri top100 CG count density", "adjust=5" ))
parkeri.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "parkeri", max.count = 150,adjust = 5 ) + labs(title= paste("parkeri top100 CG count density", "adjust=5" ))
parkeri.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "parkeri", max.count = 150,adjust = 5 ) + labs(title= paste("parkeri top100 CG count density", "adjust=5" ))

parkeri.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "parkeri", max.count = 150,adjust = 10) + labs(title= paste("parkeri top100 CG count density", "adjust=10"))
parkeri.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "parkeri", max.count = 150,adjust = 10) + labs(title= paste("parkeri top100 CG count density", "adjust=10"))
parkeri.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "parkeri", max.count = 150,adjust = 10) + labs(title= paste("parkeri top100 CG count density", "adjust=10"))
parkeri.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "parkeri", max.count = 150,adjust = 10) + labs(title= paste("parkeri top100 CG count density", "adjust=10"))
parkeri.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "parkeri", max.count = 150,adjust = 10) + labs(title= paste("parkeri top100 CG count density", "adjust=10"))

pdf("finalPlot/parkeri-top100-pattern-count.density.pdf")

  parkeri.density.CG.ad1
  parkeri.density.AT.ad1
  parkeri.density.GC.ad1
  parkeri.density.TA.ad1
  parkeri.density.AC.ad1
  
  parkeri.density.CG.ad5
  parkeri.density.AT.ad5
  parkeri.density.GC.ad5
  parkeri.density.TA.ad5
  parkeri.density.AC.ad5
  
  parkeri.density.CG.ad10
  parkeri.density.AT.ad10
  parkeri.density.GC.ad10
  parkeri.density.TA.ad10
  parkeri.density.AC.ad10

dev.off()

## danRer7
library(BSgenome.Drerio.UCSC.danRer7)
danRer7 <- BSgenome.Drerio.UCSC.danRer7

chr="chr1"
hilbert.level <- 8
chr.seq <- danRer7[[chr]]
## adjust=1
danRer7.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"CG count density", "adjust=1" ))
danRer7.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"AT count density", "adjust=1" ))
danRer7.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"GC count density", "adjust=1" ))
danRer7.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"TA count density", "adjust=1" ))
danRer7.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"AC count density", "adjust=1" ))
## adjust=5
danRer7.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"CG count density", "adjust=5" ))
danRer7.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"AT count density", "adjust=5" ))
danRer7.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"GC count density", "adjust=5" ))
danRer7.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"TA count density", "adjust=5" ))
danRer7.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"AC count density", "adjust=5" ))
## adjust=10
danRer7.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"CG count density", "adjust=10"))
danRer7.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"AT count density", "adjust=10"))
danRer7.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"GC count density", "adjust=10"))
danRer7.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"TA count density", "adjust=10"))
danRer7.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/danRer7-",chr,"-pattern-count.density.pdf",sep = ""))

  danRer7.density.CG.ad1
  danRer7.density.AT.ad1
  danRer7.density.GC.ad1
  danRer7.density.TA.ad1
  danRer7.density.AC.ad1
  
  danRer7.density.CG.ad5
  danRer7.density.AT.ad5
  danRer7.density.GC.ad5
  danRer7.density.TA.ad5
  danRer7.density.AC.ad5
  
  danRer7.density.CG.ad10
  danRer7.density.AT.ad10
  danRer7.density.GC.ad10
  danRer7.density.TA.ad10
  danRer7.density.AC.ad10

dev.off()

chr="chr2"
hilbert.level <- 8
chr.seq <- danRer7[[chr]]
## adjust=1
danRer7.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"CG count density", "adjust=1" ))
danRer7.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"AT count density", "adjust=1" ))
danRer7.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"GC count density", "adjust=1" ))
danRer7.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"TA count density", "adjust=1" ))
danRer7.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"AC count density", "adjust=1" ))
## adjust=5
danRer7.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"CG count density", "adjust=5" ))
danRer7.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"AT count density", "adjust=5" ))
danRer7.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"GC count density", "adjust=5" ))
danRer7.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"TA count density", "adjust=5" ))
danRer7.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"AC count density", "adjust=5" ))
## adjust=10
danRer7.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"CG count density", "adjust=10"))
danRer7.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"AT count density", "adjust=10"))
danRer7.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"GC count density", "adjust=10"))
danRer7.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"TA count density", "adjust=10"))
danRer7.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/danRer7-",chr,"-pattern-count.density.pdf",sep = ""))

  danRer7.density.CG.ad1
  danRer7.density.AT.ad1
  danRer7.density.GC.ad1
  danRer7.density.TA.ad1
  danRer7.density.AC.ad1
  
  danRer7.density.CG.ad5
  danRer7.density.AT.ad5
  danRer7.density.GC.ad5
  danRer7.density.TA.ad5
  danRer7.density.AC.ad5
  
  danRer7.density.CG.ad10
  danRer7.density.AT.ad10
  danRer7.density.GC.ad10
  danRer7.density.TA.ad10
  danRer7.density.AC.ad10

dev.off()

chr="chr3"
hilbert.level <- 8
chr.seq <- danRer7[[chr]]
## adjust=1
danRer7.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"CG count density", "adjust=1" ))
danRer7.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"AT count density", "adjust=1" ))
danRer7.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"GC count density", "adjust=1" ))
danRer7.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"TA count density", "adjust=1" ))
danRer7.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"AC count density", "adjust=1" ))
## adjust=5
danRer7.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"CG count density", "adjust=5" ))
danRer7.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"AT count density", "adjust=5" ))
danRer7.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"GC count density", "adjust=5" ))
danRer7.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"TA count density", "adjust=5" ))
danRer7.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"AC count density", "adjust=5" ))
## adjust=10
danRer7.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"CG count density", "adjust=10"))
danRer7.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"AT count density", "adjust=10"))
danRer7.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"GC count density", "adjust=10"))
danRer7.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"TA count density", "adjust=10"))
danRer7.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/danRer7-",chr,"-pattern-count.density.pdf",sep = ""))

  danRer7.density.CG.ad1
  danRer7.density.AT.ad1
  danRer7.density.GC.ad1
  danRer7.density.TA.ad1
  danRer7.density.AC.ad1
  
  danRer7.density.CG.ad5
  danRer7.density.AT.ad5
  danRer7.density.GC.ad5
  danRer7.density.TA.ad5
  danRer7.density.AC.ad5
  
  danRer7.density.CG.ad10
  danRer7.density.AT.ad10
  danRer7.density.GC.ad10
  danRer7.density.TA.ad10
  danRer7.density.AC.ad10

dev.off()

chr="chr4"
hilbert.level <- 8
chr.seq <- danRer7[[chr]]
## adjust=1
danRer7.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"CG count density", "adjust=1" ))
danRer7.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"AT count density", "adjust=1" ))
danRer7.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"GC count density", "adjust=1" ))
danRer7.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"TA count density", "adjust=1" ))
danRer7.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"AC count density", "adjust=1" ))
## adjust=5
danRer7.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"CG count density", "adjust=5" ))
danRer7.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"AT count density", "adjust=5" ))
danRer7.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"GC count density", "adjust=5" ))
danRer7.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"TA count density", "adjust=5" ))
danRer7.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"AC count density", "adjust=5" ))
## adjust=10
danRer7.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"CG count density", "adjust=10"))
danRer7.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"AT count density", "adjust=10"))
danRer7.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"GC count density", "adjust=10"))
danRer7.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"TA count density", "adjust=10"))
danRer7.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/danRer7-",chr,"-pattern-count.density.pdf",sep = ""))

  danRer7.density.CG.ad1
  danRer7.density.AT.ad1
  danRer7.density.GC.ad1
  danRer7.density.TA.ad1
  danRer7.density.AC.ad1
  
  danRer7.density.CG.ad5
  danRer7.density.AT.ad5
  danRer7.density.GC.ad5
  danRer7.density.TA.ad5
  danRer7.density.AC.ad5
  
  danRer7.density.CG.ad10
  danRer7.density.AT.ad10
  danRer7.density.GC.ad10
  danRer7.density.TA.ad10
  danRer7.density.AC.ad10

dev.off()

chr="chr5"
hilbert.level <- 8
chr.seq <- danRer7[[chr]]
## adjust=1
danRer7.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"CG count density", "adjust=1" ))
danRer7.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"AT count density", "adjust=1" ))
danRer7.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"GC count density", "adjust=1" ))
danRer7.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"TA count density", "adjust=1" ))
danRer7.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"AC count density", "adjust=1" ))
## adjust=5
danRer7.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"CG count density", "adjust=5" ))
danRer7.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"AT count density", "adjust=5" ))
danRer7.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"GC count density", "adjust=5" ))
danRer7.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"TA count density", "adjust=5" ))
danRer7.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"AC count density", "adjust=5" ))
## adjust=10
danRer7.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"CG count density", "adjust=10"))
danRer7.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"AT count density", "adjust=10"))
danRer7.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"GC count density", "adjust=10"))
danRer7.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"TA count density", "adjust=10"))
danRer7.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/danRer7-",chr,"-pattern-count.density.pdf",sep = ""))

  danRer7.density.CG.ad1
  danRer7.density.AT.ad1
  danRer7.density.GC.ad1
  danRer7.density.TA.ad1
  danRer7.density.AC.ad1
  
  danRer7.density.CG.ad5
  danRer7.density.AT.ad5
  danRer7.density.GC.ad5
  danRer7.density.TA.ad5
  danRer7.density.AC.ad5
  
  danRer7.density.CG.ad10
  danRer7.density.AT.ad10
  danRer7.density.GC.ad10
  danRer7.density.TA.ad10
  danRer7.density.AC.ad10

dev.off()

chr="chr19"
hilbert.level <- 7
chr.seq <- danRer7[[chr]]
## adjust=1
danRer7.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"CG count density", "adjust=1" ))
danRer7.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"AT count density", "adjust=1" ))
danRer7.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"GC count density", "adjust=1" ))
danRer7.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"TA count density", "adjust=1" ))
danRer7.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 1 ) + labs(title= paste("danRer7",chr,"AC count density", "adjust=1" ))
## adjust=5
danRer7.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"CG count density", "adjust=5" ))
danRer7.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"AT count density", "adjust=5" ))
danRer7.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"GC count density", "adjust=5" ))
danRer7.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"TA count density", "adjust=5" ))
danRer7.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 5 ) + labs(title= paste("danRer7",chr,"AC count density", "adjust=5" ))
## adjust=10
danRer7.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"CG count density", "adjust=10"))
danRer7.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"AT count density", "adjust=10"))
danRer7.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"GC count density", "adjust=10"))
danRer7.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"TA count density", "adjust=10"))
danRer7.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "danRer7", max.count = 150,adjust = 10) + labs(title= paste("danRer7",chr,"AC count density", "adjust=10"))

pdf(paste("finalPlot/danRer7-",chr,"-pattern-count.density.pdf",sep = ""))

  danRer7.density.CG.ad1
  danRer7.density.AT.ad1
  danRer7.density.GC.ad1
  danRer7.density.TA.ad1
  danRer7.density.AC.ad1
  
  danRer7.density.CG.ad5
  danRer7.density.AT.ad5
  danRer7.density.GC.ad5
  danRer7.density.TA.ad5
  danRer7.density.AC.ad5
  
  danRer7.density.CG.ad10
  danRer7.density.AT.ad10
  danRer7.density.GC.ad10
  danRer7.density.TA.ad10
  danRer7.density.AC.ad10

dev.off()

## dm6
library(BSgenome.Dmelanogaster.UCSC.dm6)
dm6 <- BSgenome.Dmelanogaster.UCSC.dm6
hilbert.level <- 9
chr.seq <- DNAString()
for(tg.seq in names(dm6)[!grepl("_", names(dm6))] ){
  chr.seq <- DNAString(c(chr.seq,dm6[[tg.seq]]))
}

dm6.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "dm6", max.count = 150,adjust = 1 ) + labs(title= paste("dm6 all CG count density", "adjust=1" ))
dm6.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "dm6", max.count = 150,adjust = 1 ) + labs(title= paste("dm6 all CG count density", "adjust=1" ))
dm6.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "dm6", max.count = 150,adjust = 1 ) + labs(title= paste("dm6 all CG count density", "adjust=1" ))
dm6.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "dm6", max.count = 150,adjust = 1 ) + labs(title= paste("dm6 all CG count density", "adjust=1" ))
dm6.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "dm6", max.count = 150,adjust = 1 ) + labs(title= paste("dm6 all CG count density", "adjust=1" ))

dm6.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "dm6", max.count = 150,adjust = 5 ) + labs(title= paste("dm6 all CG count density", "adjust=5" ))
dm6.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "dm6", max.count = 150,adjust = 5 ) + labs(title= paste("dm6 all CG count density", "adjust=5" ))
dm6.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "dm6", max.count = 150,adjust = 5 ) + labs(title= paste("dm6 all CG count density", "adjust=5" ))
dm6.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "dm6", max.count = 150,adjust = 5 ) + labs(title= paste("dm6 all CG count density", "adjust=5" ))
dm6.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "dm6", max.count = 150,adjust = 5 ) + labs(title= paste("dm6 all CG count density", "adjust=5" ))

dm6.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "dm6", max.count = 150,adjust = 10) + labs(title= paste("dm6 all CG count density", "adjust=10"))
dm6.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "dm6", max.count = 150,adjust = 10) + labs(title= paste("dm6 all CG count density", "adjust=10"))
dm6.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "dm6", max.count = 150,adjust = 10) + labs(title= paste("dm6 all CG count density", "adjust=10"))
dm6.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "dm6", max.count = 150,adjust = 10) + labs(title= paste("dm6 all CG count density", "adjust=10"))
dm6.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "dm6", max.count = 150,adjust = 10) + labs(title= paste("dm6 all CG count density", "adjust=10"))

pdf("finalPlot/dm6-all-pattern-count.density.pdf")

  dm6.density.CG.ad1
  dm6.density.AT.ad1
  dm6.density.GC.ad1
  dm6.density.TA.ad1
  dm6.density.AC.ad1
  
  dm6.density.CG.ad5
  dm6.density.AT.ad5
  dm6.density.GC.ad5
  dm6.density.TA.ad5
  dm6.density.AC.ad5
  
  dm6.density.CG.ad10
  dm6.density.AT.ad10
  dm6.density.GC.ad10
  dm6.density.TA.ad10
  dm6.density.AC.ad10

dev.off()


## ce10
library(BSgenome.Celegans.UCSC.ce10)
ce10 <- BSgenome.Celegans.UCSC.ce10
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(ce10) ){
  chr.seq <- DNAString(c(chr.seq,ce10[[tg.seq]]))
}
ce10.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "ce10", max.count = 150,adjust = 1 ) + labs(title= paste("ce10 all CG count density", "adjust=1" ))
ce10.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "ce10", max.count = 150,adjust = 1 ) + labs(title= paste("ce10 all CG count density", "adjust=1" ))
ce10.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "ce10", max.count = 150,adjust = 1 ) + labs(title= paste("ce10 all CG count density", "adjust=1" ))
ce10.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "ce10", max.count = 150,adjust = 1 ) + labs(title= paste("ce10 all CG count density", "adjust=1" ))
ce10.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "ce10", max.count = 150,adjust = 1 ) + labs(title= paste("ce10 all CG count density", "adjust=1" ))

ce10.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "ce10", max.count = 150,adjust = 5 ) + labs(title= paste("ce10 all CG count density", "adjust=5" ))
ce10.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "ce10", max.count = 150,adjust = 5 ) + labs(title= paste("ce10 all CG count density", "adjust=5" ))
ce10.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "ce10", max.count = 150,adjust = 5 ) + labs(title= paste("ce10 all CG count density", "adjust=5" ))
ce10.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "ce10", max.count = 150,adjust = 5 ) + labs(title= paste("ce10 all CG count density", "adjust=5" ))
ce10.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "ce10", max.count = 150,adjust = 5 ) + labs(title= paste("ce10 all CG count density", "adjust=5" ))

ce10.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "ce10", max.count = 150,adjust = 10) + labs(title= paste("ce10 all CG count density", "adjust=10"))
ce10.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "ce10", max.count = 150,adjust = 10) + labs(title= paste("ce10 all CG count density", "adjust=10"))
ce10.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "ce10", max.count = 150,adjust = 10) + labs(title= paste("ce10 all CG count density", "adjust=10"))
ce10.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "ce10", max.count = 150,adjust = 10) + labs(title= paste("ce10 all CG count density", "adjust=10"))
ce10.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "ce10", max.count = 150,adjust = 10) + labs(title= paste("ce10 all CG count density", "adjust=10"))

pdf("finalPlot/ce10-all-pattern-count.density.pdf")

ce10.density.CG.ad1
ce10.density.AT.ad1
ce10.density.GC.ad1
ce10.density.TA.ad1
ce10.density.AC.ad1

ce10.density.CG.ad5
ce10.density.AT.ad5
ce10.density.GC.ad5
ce10.density.TA.ad5
ce10.density.AC.ad5

ce10.density.CG.ad10
ce10.density.AT.ad10
ce10.density.GC.ad10
ce10.density.TA.ad10
ce10.density.AC.ad10

dev.off()

## tg7
library(BSgenome.Tgondii.ToxoDB.7.0)
tg7 <- Tgondii
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(tg7) ){
  chr.seq <- DNAString(c(chr.seq,tg7[[tg.seq]]))
}
tg7.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "tg7", max.count = 150,adjust = 1 ) + labs(title= paste("tg7 all CG count density", "adjust=1" ))
tg7.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "tg7", max.count = 150,adjust = 1 ) + labs(title= paste("tg7 all CG count density", "adjust=1" ))
tg7.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "tg7", max.count = 150,adjust = 1 ) + labs(title= paste("tg7 all CG count density", "adjust=1" ))
tg7.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "tg7", max.count = 150,adjust = 1 ) + labs(title= paste("tg7 all CG count density", "adjust=1" ))
tg7.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "tg7", max.count = 150,adjust = 1 ) + labs(title= paste("tg7 all CG count density", "adjust=1" ))

tg7.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "tg7", max.count = 150,adjust = 5 ) + labs(title= paste("tg7 all CG count density", "adjust=5" ))
tg7.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "tg7", max.count = 150,adjust = 5 ) + labs(title= paste("tg7 all CG count density", "adjust=5" ))
tg7.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "tg7", max.count = 150,adjust = 5 ) + labs(title= paste("tg7 all CG count density", "adjust=5" ))
tg7.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "tg7", max.count = 150,adjust = 5 ) + labs(title= paste("tg7 all CG count density", "adjust=5" ))
tg7.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "tg7", max.count = 150,adjust = 5 ) + labs(title= paste("tg7 all CG count density", "adjust=5" ))

tg7.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "tg7", max.count = 150,adjust = 10) + labs(title= paste("tg7 all CG count density", "adjust=10"))
tg7.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "tg7", max.count = 150,adjust = 10) + labs(title= paste("tg7 all CG count density", "adjust=10"))
tg7.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "tg7", max.count = 150,adjust = 10) + labs(title= paste("tg7 all CG count density", "adjust=10"))
tg7.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "tg7", max.count = 150,adjust = 10) + labs(title= paste("tg7 all CG count density", "adjust=10"))
tg7.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "tg7", max.count = 150,adjust = 10) + labs(title= paste("tg7 all CG count density", "adjust=10"))

pdf("finalPlot/tg7-all-pattern-count.density.pdf")

tg7.density.CG.ad1
tg7.density.AT.ad1
tg7.density.GC.ad1
tg7.density.TA.ad1
tg7.density.AC.ad1

tg7.density.CG.ad5
tg7.density.AT.ad5
tg7.density.GC.ad5
tg7.density.TA.ad5
tg7.density.AC.ad5

tg7.density.CG.ad10
tg7.density.AT.ad10
tg7.density.GC.ad10
tg7.density.TA.ad10
tg7.density.AC.ad10

dev.off()

## sacCer3
library(BSgenome.Scerevisiae.UCSC.sacCer3)
saccer3<- Scerevisiae
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(saccer3) ){
  chr.seq <- DNAString(c(chr.seq,saccer3[[tg.seq]]))
}
sacCer3.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "sacCer3", max.count = 150,adjust = 1 ) + labs(title= paste("sacCer3 all CG count density", "adjust=1" ))
sacCer3.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "sacCer3", max.count = 150,adjust = 1 ) + labs(title= paste("sacCer3 all CG count density", "adjust=1" ))
sacCer3.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "sacCer3", max.count = 150,adjust = 1 ) + labs(title= paste("sacCer3 all CG count density", "adjust=1" ))
sacCer3.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "sacCer3", max.count = 150,adjust = 1 ) + labs(title= paste("sacCer3 all CG count density", "adjust=1" ))
sacCer3.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "sacCer3", max.count = 150,adjust = 1 ) + labs(title= paste("sacCer3 all CG count density", "adjust=1" ))

sacCer3.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "sacCer3", max.count = 150,adjust = 5 ) + labs(title= paste("sacCer3 all CG count density", "adjust=5" ))
sacCer3.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "sacCer3", max.count = 150,adjust = 5 ) + labs(title= paste("sacCer3 all CG count density", "adjust=5" ))
sacCer3.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "sacCer3", max.count = 150,adjust = 5 ) + labs(title= paste("sacCer3 all CG count density", "adjust=5" ))
sacCer3.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "sacCer3", max.count = 150,adjust = 5 ) + labs(title= paste("sacCer3 all CG count density", "adjust=5" ))
sacCer3.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "sacCer3", max.count = 150,adjust = 5 ) + labs(title= paste("sacCer3 all CG count density", "adjust=5" ))

sacCer3.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "sacCer3", max.count = 150,adjust = 10) + labs(title= paste("sacCer3 all CG count density", "adjust=10"))
sacCer3.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "sacCer3", max.count = 150,adjust = 10) + labs(title= paste("sacCer3 all CG count density", "adjust=10"))
sacCer3.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "sacCer3", max.count = 150,adjust = 10) + labs(title= paste("sacCer3 all CG count density", "adjust=10"))
sacCer3.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "sacCer3", max.count = 150,adjust = 10) + labs(title= paste("sacCer3 all CG count density", "adjust=10"))
sacCer3.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "sacCer3", max.count = 150,adjust = 10) + labs(title= paste("sacCer3 all CG count density", "adjust=10"))

pdf("finalPlot/sacCer3-all-pattern-count.density.pdf")

sacCer3.density.CG.ad1
sacCer3.density.AT.ad1
sacCer3.density.GC.ad1
sacCer3.density.TA.ad1
sacCer3.density.AC.ad1

sacCer3.density.CG.ad5
sacCer3.density.AT.ad5
sacCer3.density.GC.ad5
sacCer3.density.TA.ad5
sacCer3.density.AC.ad5

sacCer3.density.CG.ad10
sacCer3.density.AT.ad10
sacCer3.density.GC.ad10
sacCer3.density.TA.ad10
sacCer3.density.AC.ad10

dev.off()

## ecoli
library(BSgenome.Ecoli.NCBI.20080805)
ecoli <- BSgenome.Ecoli.NCBI.20080805
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(ecoli) ){
  chr.seq <- DNAString(c(chr.seq,ecoli[[tg.seq]]))
}
ecoli.density.CG.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "ecoli", max.count = 150,adjust = 1 ) + labs(title= paste("ecoli all CG count density", "adjust=1" ))
ecoli.density.AT.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "ecoli", max.count = 150,adjust = 1 ) + labs(title= paste("ecoli all CG count density", "adjust=1" ))
ecoli.density.GC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "ecoli", max.count = 150,adjust = 1 ) + labs(title= paste("ecoli all CG count density", "adjust=1" ))
ecoli.density.TA.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "ecoli", max.count = 150,adjust = 1 ) + labs(title= paste("ecoli all CG count density", "adjust=1" ))
ecoli.density.AC.ad1  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "ecoli", max.count = 150,adjust = 1 ) + labs(title= paste("ecoli all CG count density", "adjust=1" ))

ecoli.density.CG.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "ecoli", max.count = 150,adjust = 5 ) + labs(title= paste("ecoli all CG count density", "adjust=5" ))
ecoli.density.AT.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "ecoli", max.count = 150,adjust = 5 ) + labs(title= paste("ecoli all CG count density", "adjust=5" ))
ecoli.density.GC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "ecoli", max.count = 150,adjust = 5 ) + labs(title= paste("ecoli all CG count density", "adjust=5" ))
ecoli.density.TA.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "ecoli", max.count = 150,adjust = 5 ) + labs(title= paste("ecoli all CG count density", "adjust=5" ))
ecoli.density.AC.ad5  <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "ecoli", max.count = 150,adjust = 5 ) + labs(title= paste("ecoli all CG count density", "adjust=5" ))

ecoli.density.CG.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "CG",  species = "ecoli", max.count = 150,adjust = 10) + labs(title= paste("ecoli all CG count density", "adjust=10"))
ecoli.density.AT.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AT",  species = "ecoli", max.count = 150,adjust = 10) + labs(title= paste("ecoli all CG count density", "adjust=10"))
ecoli.density.GC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "GC",  species = "ecoli", max.count = 150,adjust = 10) + labs(title= paste("ecoli all CG count density", "adjust=10"))
ecoli.density.TA.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "TA",  species = "ecoli", max.count = 150,adjust = 10) + labs(title= paste("ecoli all CG count density", "adjust=10"))
ecoli.density.AC.ad10 <- patternDensityPlot(chr.seq,level=hilbert.level, pattern = "AC",  species = "ecoli", max.count = 150,adjust = 10) + labs(title= paste("ecoli all CG count density", "adjust=10"))

pdf("finalPlot/ecoli-all-pattern-count.density.pdf")

ecoli.density.CG.ad1
ecoli.density.AT.ad1
ecoli.density.GC.ad1
ecoli.density.TA.ad1
ecoli.density.AC.ad1

ecoli.density.CG.ad5
ecoli.density.AT.ad5
ecoli.density.GC.ad5
ecoli.density.TA.ad5
ecoli.density.AC.ad5

ecoli.density.CG.ad10
ecoli.density.AT.ad10
ecoli.density.GC.ad10
ecoli.density.TA.ad10
ecoli.density.AC.ad10

dev.off()


# pdf("finalPlot/CG.density.pdf",height = 3, width = 5)
# hg38.hist.CG
# mm10.hist.CG
# taeGut1.hist.CG
# parkeri.hist.CG
# danRer7.hist.CG
# dm6.hist.CG
# ce10.hist.CG
# tg7.hist.AT
# sacCer3.hist.CG
# ecoli.hist.CG
# dev.off()
# 
# pdf("finalPlot/AT.desity.pdf",height = 3, width = 5)
# hg38.hist.AT
# mm10.hist.AT
# taeGut1.hist.AT
# parkeri.hist.AT
# danRer7.hist.AT
# dm6.hist.AT
# ce10.hist.AT
# tg7.hist.AT
# sacCer3.hist.AT
# ecoli.hist.AT
# dev.off()

## plot wig
##hg38
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38

chr="chr1"
chr.seq <- hg38[[chr]]
# wigPlot(chr.seq,pattern = "CG",wig.file = paste("finalPlot/hg38-",chr,"-CG-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
# wigPlot(chr.seq,pattern = "AT",wig.file = paste("finalPlot/hg38-",chr,"-AT-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "GC",wig.file = paste("finalPlot/hg38-",chr,"-GC-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "TA",wig.file = paste("finalPlot/hg38-",chr,"-TA-win2500-ps.wig",sep=""),win = 1250,chrom = chr)

chr="chr2"
chr.seq <- hg38[[chr]]
wigPlot(chr.seq,pattern = "CG",wig.file = paste("finalPlot/hg38-",chr,"-CG-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "AT",wig.file = paste("finalPlot/hg38-",chr,"-AT-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "GC",wig.file = paste("finalPlot/hg38-",chr,"-GC-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "TA",wig.file = paste("finalPlot/hg38-",chr,"-TA-win2500-ps.wig",sep=""),win = 1250,chrom = chr)

chr="chr3"
chr.seq <- hg38[[chr]]
wigPlot(chr.seq,pattern = "CG",wig.file = paste("finalPlot/hg38-",chr,"-CG-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "AT",wig.file = paste("finalPlot/hg38-",chr,"-AT-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "GC",wig.file = paste("finalPlot/hg38-",chr,"-GC-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "TA",wig.file = paste("finalPlot/hg38-",chr,"-TA-win2500-ps.wig",sep=""),win = 1250,chrom = chr)

chr="chr3"
chr.seq <- hg38[[chr]]
wigPlot(chr.seq,pattern = "CG",wig.file = paste("finalPlot/hg38-",chr,"-CG-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "AT",wig.file = paste("finalPlot/hg38-",chr,"-AT-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "GC",wig.file = paste("finalPlot/hg38-",chr,"-GC-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "TA",wig.file = paste("finalPlot/hg38-",chr,"-TA-win2500-ps.wig",sep=""),win = 1250,chrom = chr)



##mm10
library(BSgenome.Mmusculus.UCSC.mm10)
mm10 <- BSgenome.Mmusculus.UCSC.mm10

chr="chr1"
chr.seq <- mm10[[chr]]
# wigPlot(chr.seq,pattern = "CG",wig.file = paste("finalPlot/mm10-",chr,"-CG-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
# wigPlot(chr.seq,pattern = "AT",wig.file = paste("finalPlot/mm10-",chr,"-AT-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "GC",wig.file = paste("finalPlot/mm10-",chr,"-GC-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "TA",wig.file = paste("finalPlot/mm10-",chr,"-TA-win2500-ps.wig",sep=""),win = 1250,chrom = chr)

chr="chr2"
chr.seq <- mm10[[chr]]
wigPlot(chr.seq,pattern = "CG",wig.file = paste("finalPlot/mm10-",chr,"-CG-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "AT",wig.file = paste("finalPlot/mm10-",chr,"-AT-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "GC",wig.file = paste("finalPlot/mm10-",chr,"-GC-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "TA",wig.file = paste("finalPlot/mm10-",chr,"-TA-win2500-ps.wig",sep=""),win = 1250,chrom = chr)

chr="chr3"
chr.seq <- mm10[[chr]]
wigPlot(chr.seq,pattern = "CG",wig.file = paste("finalPlot/mm10-",chr,"-CG-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "AT",wig.file = paste("finalPlot/mm10-",chr,"-AT-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "GC",wig.file = paste("finalPlot/mm10-",chr,"-GC-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "TA",wig.file = paste("finalPlot/mm10-",chr,"-TA-win2500-ps.wig",sep=""),win = 1250,chrom = chr)

chr="chr3"
chr.seq <- mm10[[chr]]
wigPlot(chr.seq,pattern = "CG",wig.file = paste("finalPlot/mm10-",chr,"-CG-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "AT",wig.file = paste("finalPlot/mm10-",chr,"-AT-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "GC",wig.file = paste("finalPlot/mm10-",chr,"-GC-win2500-ps.wig",sep=""),win = 1250,chrom = chr)
wigPlot(chr.seq,pattern = "TA",wig.file = paste("finalPlot/mm10-",chr,"-TA-win2500-ps.wig",sep=""),win = 1250,chrom = chr)


## ce10
library(BSgenome.Celegans.UCSC.ce10)
ce10 <- BSgenome.Celegans.UCSC.ce10
chr.seq <- DNAString()
for(tg.seq in names(ce10) ){
  chr.seq <- DNAString(c(chr.seq,ce10[[tg.seq]]))
}
# wigPlot(chr.seq,pattern = "CG",wig.file = "finalPlot/ce10-chr1-CG-win2500-ps.wig",win = 1250,chrom = "chr1")
# wigPlot(chr.seq,pattern = "AT",wig.file = "finalPlot/ce10-chr1-AT-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "GC",wig.file = "finalPlot/ce10-chr1-GC-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "TA",wig.file = "finalPlot/ce10-chr1-TA-win2500-ps.wig",win = 1250,chrom = "chr1")

## danRer7
library(BSgenome.Drerio.UCSC.danRer7)
danRer7 <- BSgenome.Drerio.UCSC.danRer7
chr.seq <- danRer7[[1]]
# wigPlot(chr.seq,pattern = "CG",wig.file = "finalPlot/danRer7-chr1-CG-win2500-ps.wig",win = 1250,chrom = "chr1")
# wigPlot(chr.seq,pattern = "AT",wig.file = "finalPlot/danRer7-chr1-AT-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "GC",wig.file = "finalPlot/danRer7-chr1-GC-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "TA",wig.file = "finalPlot/danRer7-chr1-TA-win2500-ps.wig",win = 1250,chrom = "chr1")

## dm6
library(BSgenome.Dmelanogaster.UCSC.dm6)
dm6 <- BSgenome.Dmelanogaster.UCSC.dm6
chr.seq <- DNAString()
for(tg.seq in names(dm6)[!grepl("_", names(dm6))] ){
  chr.seq <- DNAString(c(chr.seq,dm6[[tg.seq]]))
}
# wigPlot(chr.seq,pattern = "CG",wig.file = "finalPlot/dm6-chr1-CG-win2500-ps.wig",win = 1250,chrom = "chr1")
# wigPlot(chr.seq,pattern = "AT",wig.file = "finalPlot/dm6-chr1-AT-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "GC",wig.file = "finalPlot/dm6-chr1-GC-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "GC",wig.file = "finalPlot/dm6-chr1-TA-win2500-ps.wig",win = 1250,chrom = "chr1")

## ecoli
library(BSgenome.Ecoli.NCBI.20080805)
ecoli <- BSgenome.Ecoli.NCBI.20080805
chr.seq <- DNAString()
for(tg.seq in names(ecoli) ){
  chr.seq <- DNAString(c(chr.seq,ecoli[[tg.seq]]))
}
# wigPlot(chr.seq,pattern = "CG",wig.file = "finalPlot/ecoli-chr1-CG-win2500-ps.wig",win = 1250,chrom = "chr1")
# wigPlot(chr.seq,pattern = "AT",wig.file = "finalPlot/ecoli-chr1-AT-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "GC",wig.file = "finalPlot/ecoli-chr1-GC-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "TA",wig.file = "finalPlot/ecoli-chr1-TA-win2500-ps.wig",win = 1250,chrom = "chr1")

## parkeri
library("Biostrings")
parkeri.path <- "~/projects/CpGden/parkeri.fasta"
parkeri <- readDNAStringSet(parkeri.path, format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=F)
hilbert.level <- 9
pkr <- DNAString()
for(tg.seq in c(1:100) ){
  # pkr <- DNAString(c(pkr,parkeri[tg.seq]))
  pkr <- DNAStringSet(paste(as.character(pkr),as.character(parkeri[tg.seq]),sep=""))
}
chr.seq <- pkr[[1]]
# wigPlot(chr.seq,pattern = "CG",wig.file = "finalPlot/parkeri-chr1-CG-win2500-ps.wig",win = 1250,chrom = "chr1")
# wigPlot(chr.seq,pattern = "AT",wig.file = "finalPlot/parkeri-chr1-AT-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "GC",wig.file = "finalPlot/parkeri-chr1-GC-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "TA",wig.file = "finalPlot/parkeri-chr1-TA-win2500-ps.wig",win = 1250,chrom = "chr1")

## sacCer3
library(BSgenome.Scerevisiae.UCSC.sacCer3)
saccer3<- Scerevisiae
hilbert.level <- 7
chr.seq <- DNAString()
for(tg.seq in names(saccer3) ){
  chr.seq <- DNAString(c(chr.seq,saccer3[[tg.seq]]))
}
# wigPlot(chr.seq,pattern = "CG",wig.file = "finalPlot/sacCer3-chr1-CG-win2500-ps.wig",win = 1250,chrom = "chr1")
# wigPlot(chr.seq,pattern = "AT",wig.file = "finalPlot/sacCer3-chr1-AT-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "GC",wig.file = "finalPlot/sacCer3-chr1-GC-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "TA",wig.file = "finalPlot/sacCer3-chr1-TA-win2500-ps.wig",win = 1250,chrom = "chr1")

## taeGut1
library(BSgenome.Tguttata.UCSC.taeGut1)
taeGut1 <- BSgenome.Tguttata.UCSC.taeGut1
chr.seq <- taeGut1[[1]]
# wigPlot(chr.seq,pattern = "CG",wig.file = "finalPlot/taeGut1-chr1-CG-win2500-ps.wig",win = 1250,chrom = "chr1")
# wigPlot(chr.seq,pattern = "AT",wig.file = "finalPlot/taeGut1-chr1-AT-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "GC",wig.file = "finalPlot/taeGut1-chr1-GC-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "TA",wig.file = "finalPlot/taeGut1-chr1-TA-win2500-ps.wig",win = 1250,chrom = "chr1")

## tg7
library(BSgenome.Tgondii.ToxoDB.7.0)
tg7 <- Tgondii
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(tg7) ){
  chr.seq <- DNAString(c(chr.seq,tg7[[tg.seq]]))
}
# wigPlot(chr.seq,pattern = "CG",wig.file = "finalPlot/tg7-chr1-CG-win2500-ps.wig",win = 1250,chrom = "chr1")
# wigPlot(chr.seq,pattern = "AT",wig.file = "finalPlot/tg7-chr1-AT-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "GC",wig.file = "finalPlot/tg7-chr1-GC-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "TA",wig.file = "finalPlot/tg7-chr1-TA-win2500-ps.wig",win = 1250,chrom = "chr1")




## boxplot
##hg38
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38
hilbert.level <- 9
chr.seq <- hg38[[1]]
hg38.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
hg38.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
hg38.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
hg38.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
hg38.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")

##mm10
library(BSgenome.Mmusculus.UCSC.mm10)
mm10 <- BSgenome.Mmusculus.UCSC.mm10
hilbert.level <- 9
chr.seq <- mm10[[1]]
mm10.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
mm10.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
mm10.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
mm10.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
mm10.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")


## ce10
library(BSgenome.Celegans.UCSC.ce10)
ce10 <- BSgenome.Celegans.UCSC.ce10
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(ce10) ){
  chr.seq <- DNAString(c(chr.seq,ce10[[tg.seq]]))
}
ce10.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
ce10.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
ce10.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
ce10.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
ce10.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")


## danRer7
library(BSgenome.Drerio.UCSC.danRer7)
danRer7 <- BSgenome.Drerio.UCSC.danRer7
chr.seq <- danRer7[[1]]
hilbert.level <- 8
danRer7.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
danRer7.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
danRer7.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
danRer7.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
danRer7.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")


## dm6
library(BSgenome.Dmelanogaster.UCSC.dm6)
dm6 <- BSgenome.Dmelanogaster.UCSC.dm6
hilbert.level <- 9
chr.seq <- DNAString()
for(tg.seq in names(dm6)[!grepl("_", names(dm6))] ){
  chr.seq <- DNAString(c(chr.seq,dm6[[tg.seq]]))
}
dm6.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
dm6.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
dm6.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
dm6.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
dm6.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")

## ecoli
library(BSgenome.Ecoli.NCBI.20080805)
ecoli <- BSgenome.Ecoli.NCBI.20080805
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(ecoli) ){
  chr.seq <- DNAString(c(chr.seq,ecoli[[tg.seq]]))
}
ecoli.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
ecoli.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
ecoli.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
ecoli.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
ecoli.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")


## parkeri
library("Biostrings")
parkeri.path <- "~/projects/CpGden/parkeri.fasta"
parkeri <- readDNAStringSet(parkeri.path, format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=F)
hilbert.level <- 9
pkr <- DNAString()
for(tg.seq in c(1:100) ){
  # pkr <- DNAString(c(pkr,parkeri[tg.seq]))
  pkr <- DNAStringSet(paste(as.character(pkr),as.character(parkeri[tg.seq]),sep=""))
}
chr.seq <- pkr[[1]]
parkeri.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
parkeri.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
parkeri.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
parkeri.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
parkeri.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")


## sacCer3
library(BSgenome.Scerevisiae.UCSC.sacCer3)
saccer3<- Scerevisiae
hilbert.level <- 7
chr.seq <- DNAString()
for(tg.seq in names(saccer3) ){
  chr.seq <- DNAString(c(chr.seq,saccer3[[tg.seq]]))
}
sacCer3.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
sacCer3.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
sacCer3.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
sacCer3.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
sacCer3.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")

## taeGut1
library(BSgenome.Tguttata.UCSC.taeGut1)
taeGut1 <- BSgenome.Tguttata.UCSC.taeGut1
hilbert.level <- 8
chr.seq <- taeGut1[[1]]
taeGut1.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
taeGut1.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
taeGut1.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
taeGut1.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
taeGut1.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")

## tg7
library(BSgenome.Tgondii.ToxoDB.7.0)
tg7 <- Tgondii
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(tg7) ){
  chr.seq <- DNAString(c(chr.seq,tg7[[tg.seq]]))
}
tg7.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
tg7.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
tg7.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
tg7.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
tg7.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")


CG.specie <- c(rep("hg38",length(hg38.CG.count)),
               rep("mm10",length(mm10.CG.count)),
               rep("taeGut1",length(taeGut1.CG.count)),
               rep("parkeri",length(parkeri.CG.count)),
               rep("danRer7",length(danRer7.CG.count)),
               rep("dm6",length(dm6.CG.count)),
               rep("ce10",length(ce10.CG.count)),
               rep("tg7",length(tg7.CG.count)),
               rep("sacCer3",length(sacCer3.CG.count)),
               rep("ecoli",length(ecoli.CG.count)))

CG.count <- c(hg38.CG.count,mm10.CG.count,taeGut1.CG.count,parkeri.CG.count,danRer7.CG.count,dm6.CG.count,ce10.CG.count,tg7.CG.count,sacCer3.CG.count,ecoli.CG.count) 

AT.specie <- c(rep("hg38",length(hg38.AT.count)),
               rep("mm10",length(mm10.AT.count)),
               rep("taeGut1",length(taeGut1.AT.count)),
               rep("parkeri",length(parkeri.AT.count)),
               rep("danRer7",length(danRer7.AT.count)),
               rep("dm6",length(dm6.AT.count)),
               rep("ce10",length(ce10.AT.count)),
               rep("tg7",length(tg7.AT.count)),
               rep("sacCer3",length(sacCer3.AT.count)),
               rep("ecoli",length(ecoli.AT.count)))

AT.count <- c(hg38.AT.count,mm10.AT.count,taeGut1.AT.count,parkeri.AT.count,danRer7.AT.count,dm6.AT.count,ce10.AT.count,tg7.AT.count,sacCer3.AT.count,ecoli.AT.count) 

GC.specie <- c(rep("hg38",length(hg38.GC.count)),
               rep("mm10",length(mm10.GC.count)),
               rep("taeGut1",length(taeGut1.GC.count)),
               rep("parkeri",length(parkeri.GC.count)),
               rep("danRer7",length(danRer7.GC.count)),
               rep("dm6",length(dm6.GC.count)),
               rep("ce10",length(ce10.GC.count)),
               rep("tg7",length(tg7.GC.count)),
               rep("sacCer3",length(sacCer3.GC.count)),
               rep("ecoli",length(ecoli.GC.count)))

GC.count <- c(hg38.GC.count,mm10.GC.count,taeGut1.GC.count,parkeri.GC.count,danRer7.GC.count,dm6.GC.count,ce10.GC.count,tg7.GC.count,sacCer3.GC.count,ecoli.GC.count) 

TA.specie <- c(rep("hg38",length(hg38.TA.count)),
               rep("mm10",length(mm10.TA.count)),
               rep("taeGut1",length(taeGut1.TA.count)),
               rep("parkeri",length(parkeri.TA.count)),
               rep("danRer7",length(danRer7.TA.count)),
               rep("dm6",length(dm6.TA.count)),
               rep("ce10",length(ce10.TA.count)),
               rep("tg7",length(tg7.TA.count)),
               rep("sacCer3",length(sacCer3.TA.count)),
               rep("ecoli",length(ecoli.TA.count)))

TA.count <- c(hg38.TA.count,mm10.TA.count,taeGut1.TA.count,parkeri.TA.count,danRer7.TA.count,dm6.TA.count,ce10.TA.count,tg7.TA.count,sacCer3.TA.count,ecoli.TA.count) 

AC.specie <- c(rep("hg38",length(hg38.AC.count)),
               rep("mm10",length(mm10.AC.count)),
               rep("taeGut1",length(taeGut1.AC.count)),
               rep("parkeri",length(parkeri.AC.count)),
               rep("danRer7",length(danRer7.AC.count)),
               rep("dm6",length(dm6.AC.count)),
               rep("ce10",length(ce10.AC.count)),
               rep("tg7",length(tg7.AC.count)),
               rep("sacCer3",length(sacCer3.AC.count)),
               rep("ecoli",length(ecoli.AC.count)))

AC.count <- c(hg38.AC.count,mm10.AC.count,taeGut1.AC.count,parkeri.AC.count,danRer7.AC.count,dm6.AC.count,ce10.AC.count,tg7.AC.count,sacCer3.AC.count,ecoli.AC.count) 



bin.CG.count <- data.frame(specie=CG.specie,count=CG.count)
bin.AT.count <- data.frame(specie=AT.specie,count=AT.count)
bin.GC.count <- data.frame(specie=GC.specie,count=GC.count)
bin.TA.count <- data.frame(specie=TA.specie,count=TA.count)
bin.AC.count <- data.frame(specie=AC.specie,count=AC.count)

bin.CG.count$specie <- factor(bin.CG.count$specie,levels = c("hg38","mm10","taeGut1","parkeri","danRer7","dm6","ce10","tg7","sacCer3","ecoli"))
bin.AT.count$specie <- factor(bin.AT.count$specie,levels = c("hg38","mm10","taeGut1","parkeri","danRer7","dm6","ce10","tg7","sacCer3","ecoli"))
bin.GC.count$specie <- factor(bin.GC.count$specie,levels = c("hg38","mm10","taeGut1","parkeri","danRer7","dm6","ce10","tg7","sacCer3","ecoli"))
bin.TA.count$specie <- factor(bin.TA.count$specie,levels = c("hg38","mm10","taeGut1","parkeri","danRer7","dm6","ce10","tg7","sacCer3","ecoli"))
bin.AC.count$specie <- factor(bin.AC.count$specie,levels = c("hg38","mm10","taeGut1","parkeri","danRer7","dm6","ce10","tg7","sacCer3","ecoli"))


pdf(file = "finalPlot/all.specie.bin.pattern.count.boxplot.pdf")
ggplot(bin.CG.count,aes(x = specie, y=count)) + geom_boxplot() + labs(title="CG count / 1kbp")
ggplot(bin.AT.count,aes(x = specie, y=count)) + geom_boxplot() + labs(title="AT count / 1kbp")
ggplot(bin.GC.count,aes(x = specie, y=count)) + geom_boxplot() + labs(title="GC count / 1kbp")
ggplot(bin.TA.count,aes(x = specie, y=count)) + geom_boxplot() + labs(title="TA count / 1kbp")
ggplot(bin.AC.count,aes(x = specie, y=count)) + geom_boxplot() + labs(title="AC count / 1kbp")
dev.off()




### Hilbert

hilbertPlot <- function(chr.seq,target.pattern = "CG",hilbert.level=9,specie="hg38",chr="chr1",pdf,legbreak=c(30,20,10,0),legcolor= c("#EE0000", "#FF4500","#EE8262","#EED5B7")){
  pt <- matchPattern(target.pattern,chr.seq)
  pt.pos <- start(pt)
  hbc = HilbertCurve(1, length(chr.seq), level = hilbert.level, mode = "pixel")
  chr.brks <- hbc@BINS@start
  pt.nu <- table(cut(pt.pos,c(chr.brks,length(chr.seq))))
  pt.bin.nu <- as.integer(pt.nu)/hbc@BINS@width*1000
  pt.bin.start <- chr.brks
  pt.bin.end <- c(pt.bin.start[-1],length(chr.seq))
  pt.bin <- IRanges(pt.bin.start, pt.bin.end)
  col_fun = colorRamp2(legbreak, legcolor)
  cm= ColorMapping(levels = legbreak ,colors =legcolor)
  legend = color_mapping_legend(cm, plot = F, title =  paste(target.pattern,"count"))
  pdf(pdf)
  hc = HilbertCurve(1, length(chr.seq), level = hilbert.level, mode = "pixel",title = paste(target.pattern,"Count distribution of",specie,chr),legend = legend)
  hc_layer(hc,pt.bin,col=col_fun(pt.bin.nu))
  dev.off()
}


##hg38
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38
specie="hg38"
legcolor= c("#EE0000", "#FF4500","#EE8262","#EED5B7")

for(chr in seqnames(hg38)[c(1:5,19)]){
  hilbert.level <- 9
  chr.seq <- hg38[[chr]]
  target.pattern <- "CG"
  legbreak=c(30,20,10,0)
  pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
  hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
              specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)
  }

for(chr in seqnames(hg38)[c(1:5,19)]){
  hilbert.level <- 9
  chr.seq <- hg38[[chr]]
  target.pattern <- "AT"
  legbreak=c(100,75,50,0)
  pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
  hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)
}

## mm10

library(BSgenome.Mmusculus.UCSC.mm10)
mm10 <- BSgenome.Mmusculus.UCSC.mm10

specie="mm10"
legcolor= c("#EE0000", "#FF4500","#EE8262","#EED5B7")

for(chr in seqnames(mm10)[c(1:5,19)]){
  hilbert.level <- 9
  chr.seq <- mm10[[chr]]
  target.pattern <- "CG"
  legbreak=c(30,20,10,0)
  pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
  hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
              specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)
}

for(chr in seqnames(mm10)[c(1:5,19)]){
  hilbert.level <- 9
  chr.seq <- mm10[[chr]]
  target.pattern <- "AT"
  legbreak=c(100,75,50,0)
  pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
  hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
              specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)
}


## taeGut1
library(BSgenome.Tguttata.UCSC.taeGut1)
taeGut1 <- BSgenome.Tguttata.UCSC.taeGut1

specie="taeGut1"
legcolor= c("#EE0000", "#FF4500","#EE8262","#EED5B7")

for(chr in seqnames(taeGut1)[c(1,4,5,6,8)]){
  hilbert.level <- 9
  chr.seq <- taeGut1[[chr]]
  target.pattern <- "CG"
  legbreak=c(30,20,10,0)
  pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
  hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
              specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)
}

for(chr in seqnames(taeGut1)[c(1,4,5,6,8)]){
  hilbert.level <- 9
  chr.seq <- taeGut1[[chr]]
  target.pattern <- "AT"
  legbreak=c(100,75,50,0)
  pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
  hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
              specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)
}

## parkeri
library("Biostrings")
parkeri.path <- "~/projects/CpGden/parkeri.fasta"
parkeri <- readDNAStringSet(parkeri.path, format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=F)
hilbert.level <- 9
pkr <- DNAString()
for(tg.seq in c(1:100) ){
  # pkr <- DNAString(c(pkr,parkeri[tg.seq]))
  pkr <- DNAStringSet(paste(as.character(pkr),as.character(parkeri[tg.seq]),sep=""))
}
chr.seq <- pkr[[1]]
specie="parkeri"
target.pattern <- "CG"
legbreak=c(30,20,10,0)
pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)

target.pattern <- "AT"
legbreak=c(100,75,50,0)
pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)


## danRer7
library(BSgenome.Drerio.UCSC.danRer7)
danRer7 <- BSgenome.Drerio.UCSC.danRer7

hilbert.level <- 8
specie="danRer7"
legcolor= c("#EE0000", "#FF4500","#EE8262","#EED5B7")

for(chr in seqnames(danRer7)[c(1:5,19)]){
  chr.seq <- danRer7[[chr]]
  target.pattern <- "CG"
  legbreak=c(30,20,10,0)
  pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
  hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
              specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)
}

for(chr in seqnames(danRer7)[c(1:5,19)]){
  chr.seq <- danRer7[[chr]]
  target.pattern <- "AT"
  legbreak=c(100,75,50,0)
  pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
  hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
              specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)
}

## dm6
library(BSgenome.Dmelanogaster.UCSC.dm6)
dm6 <- BSgenome.Dmelanogaster.UCSC.dm6
hilbert.level <- 9
chr.seq <- DNAString()
for(tg.seq in names(dm6)[!grepl("_", names(dm6))] ){
  chr.seq <- DNAString(c(chr.seq,dm6[[tg.seq]]))
}

specie="dm6"
target.pattern <- "CG"
legbreak=c(30,20,10,0)
pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)

target.pattern <- "AT"
legbreak=c(100,75,50,0)
pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)



## ce10
library(BSgenome.Celegans.UCSC.ce10)
ce10 <- BSgenome.Celegans.UCSC.ce10
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(ce10) ){
  chr.seq <- DNAString(c(chr.seq,ce10[[tg.seq]]))
}
specie="ce10"
target.pattern <- "CG"
legbreak=c(30,20,10,0)
pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)

target.pattern <- "AT"
legbreak=c(100,75,50,0)
pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)


## tg7
library(BSgenome.Tgondii.ToxoDB.7.0)
tg7 <- Tgondii
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(tg7) ){
  chr.seq <- DNAString(c(chr.seq,tg7[[tg.seq]]))
}

specie="tg7"
target.pattern <- "CG"
legbreak=c(30,20,10,0)
pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)

target.pattern <- "AT"
legbreak=c(100,75,50,0)
pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)


## sacCer3
library(BSgenome.Scerevisiae.UCSC.sacCer3)
saccer3<- Scerevisiae
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(saccer3) ){
  chr.seq <- DNAString(c(chr.seq,saccer3[[tg.seq]]))
}

specie="sacCer3"
target.pattern <- "CG"
legbreak=c(30,20,10,0)
pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)

target.pattern <- "AT"
legbreak=c(100,75,50,0)
pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)


## ecoli
library(BSgenome.Ecoli.NCBI.20080805)
ecoli <- BSgenome.Ecoli.NCBI.20080805
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(ecoli) ){
  chr.seq <- DNAString(c(chr.seq,ecoli[[tg.seq]]))
}

specie="ecoli"
target.pattern <- "CG"
legbreak=c(30,20,10,0)
pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)

target.pattern <- "AT"
legbreak=c(100,75,50,0)
pdf=paste("finalPlot/hilbertPlot/",specie,"-",chr,"-",target.pattern,"-hilbertPlot.pdf",sep="")
hilbertPlot(chr.seq,target.pattern = target.pattern ,hilbert.level = hilbert.level,
            specie = specie,chr = chr,pdf = pdf,legbreak = legbreak ,legcolor = legcolor)



## distance plot
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38

chr="chr1"
chr.seq <- hg38[[chr]]
cg.pos <- matchPattern("CG",chr.seq)
cg.pos <- start(cg.pos)
cg.shif <- data.frame(pos=cg.pos, 
                      sh2 = c(cg.pos[2:length(cg.pos)],tail(cg.pos,1)),
                      sh5 = c(cg.pos[5:length(cg.pos)],tail(cg.pos,4)),
                      sh10 = c(cg.pos[10:length(cg.pos)],tail(cg.pos,9)),
                      sh15 = c(cg.pos[15:length(cg.pos)],tail(cg.pos,14)),
                      sh20 = c(cg.pos[20:length(cg.pos)],tail(cg.pos,19)),
                      sh25 = c(cg.pos[25:length(cg.pos)],tail(cg.pos,24)),
                      sh30 = c(cg.pos[30:length(cg.pos)],tail(cg.pos,29)),
                      sh35 = c(cg.pos[35:length(cg.pos)],tail(cg.pos,34)),
                      sh40 = c(cg.pos[40:length(cg.pos)],tail(cg.pos,39)),
                      sh45 = c(cg.pos[45:length(cg.pos)],tail(cg.pos,44)),
                      sh50 = c(cg.pos[50:length(cg.pos)],tail(cg.pos,49)))

cg.shif$dis2 <- with(cg.shif,sh2-pos)
cg.shif$dis5 <- with(cg.shif,sh5-pos)
cg.shif$dis10 <- with(cg.shif,sh10-pos)
cg.shif$dis15 <- with(cg.shif,sh15-pos)
cg.shif$dis20 <- with(cg.shif,sh20-pos)
cg.shif$dis25 <- with(cg.shif,sh25-pos)
cg.shif$dis30 <- with(cg.shif,sh30-pos)
cg.shif$dis35 <- with(cg.shif,sh35-pos)
cg.shif$dis40 <- with(cg.shif,sh40-pos)
cg.shif$dis45 <- with(cg.shif,sh45-pos)
cg.shif$dis50 <- with(cg.shif,sh50-pos)

CG.info <- data.table(cg.id=c(1:length(cg.pos)),pos=cg.pos,
                      dis2=cg.shif$dis2,
                      dis5=cg.shif$dis5,
                      dis10=cg.shif$dis10,
                      dis15=cg.shif$dis15,
                      dis20=cg.shif$dis20,
                      dis25=cg.shif$dis25,
                      dis30=cg.shif$dis30,
                      dis35=cg.shif$dis35,
                      dis40=cg.shif$dis40,
                      dis45=cg.shif$dis45,
                      dis50=cg.shif$dis50)



CGI <- fread("~/projects/CpGden/CpGIsland.hg38.txt")
setnames(CGI,c("chrom","start","end"))
CGI <- CGI[chrom=="chr1"]
setkey(CGI,start)
CGI[,CGI.id := c(1:nrow(CGI))]
CGI <- CGI[,.(pos=start:end),by= CGI.id]
CGI.distance <- merge(CGI,CG.info,by="pos")

pdf(file = "finalPlot/distancePlot/hg38-CGI-region-cg-distance.pdf")
ggplot(CGI.distance,aes(dis2))+geom_histogram(binwidth=1)+xlim(0,50) + labs(title="hg38 CGI region cg distance shift number = 2")
ggplot(CGI.distance)+geom_density(aes(x=dis2),fill="RoyalBlue")+xlim(0,100) + labs(title="hg38 CGI region cg distance shift number = 2")

ggplot(CGI.distance,aes(dis5))+geom_histogram(binwidth=1)+xlim(0,100) + labs(title="hg38 CGI region cg distance shift number = 5")
ggplot(CGI.distance)+geom_density(aes(x=dis5),fill="RoyalBlue")+xlim(0,100) + labs(title="hg38 CGI region cg distance shift number = 5")
ggplot(CGI.distance,aes(dis10))+geom_histogram(binwidth=1)+xlim(0,200) + labs(title="hg38 CGI region cg distance shift number = 10")
ggplot(CGI.distance)+geom_density(aes(x=dis10),fill="RoyalBlue")+xlim(0,200) + labs(title="hg38 CGI region cg distance shift number = 10")
ggplot(CGI.distance,aes(dis15))+geom_histogram(binwidth=1)+xlim(0,400) + labs(title="hg38 CGI region cg distance shift number = 15")
ggplot(CGI.distance)+geom_density(aes(x=dis15),fill="RoyalBlue")+xlim(0,400) + labs(title="hg38 CGI region cg distance shift number = 15")
ggplot(CGI.distance,aes(dis20))+geom_histogram(binwidth=1)+xlim(0,500) + labs(title="hg38 CGI region cg distance shift number = 20")
ggplot(CGI.distance)+geom_density(aes(x=dis20),fill="RoyalBlue")+xlim(0,500) + labs(title="hg38 CGI region cg distance shift number = 20")
ggplot(CGI.distance,aes(dis25))+geom_histogram(binwidth=1)+xlim(0,600) + labs(title="hg38 CGI region cg distance shift number = 25") 
ggplot(CGI.distance)+geom_density(aes(x=dis25),fill="RoyalBlue")+xlim(0,600) + labs(title="hg38 CGI region cg distance shift number = 25")
ggplot(CGI.distance,aes(dis30))+geom_histogram(binwidth=1)+xlim(0,700) + labs(title="hg38 CGI region cg distance shift number = 30")
ggplot(CGI.distance)+geom_density(aes(x=dis30),fill="RoyalBlue")+xlim(0,700) + labs(title="hg38 CGI region cg distance shift number = 30")
ggplot(CGI.distance,aes(dis35))+geom_histogram(binwidth=1)+xlim(0,1000) + labs(title="hg38 CGI region cg distance shift number = 35")
ggplot(CGI.distance)+geom_density(aes(x=dis35),fill="RoyalBlue")+xlim(0,1000) + labs(title="hg38 CGI region cg distance shift number = 35")
ggplot(CGI.distance,aes(dis40))+geom_histogram(binwidth=1)+xlim(0,1000) + labs(title="hg38 CGI region cg distance shift number = 40")
ggplot(CGI.distance)+geom_density(aes(x=dis40),fill="RoyalBlue")+xlim(0,1000) + labs(title="hg38 CGI region cg distance shift number = 40")
ggplot(CGI.distance,aes(dis45))+geom_histogram(binwidth=1)+xlim(0,1200) + labs(title="hg38 CGI region cg distance shift number = 45")
ggplot(CGI.distance)+geom_density(aes(x=dis45),fill="RoyalBlue")+xlim(0,1200) + labs(title="hg38 CGI region cg distance shift number = 45")
ggplot(CGI.distance,aes(dis50))+geom_histogram(binwidth=1)+xlim(0,3000) + labs(title="hg38 CGI region cg distance shift number = 50")
ggplot(CGI.distance)+geom_density(aes(x=dis50),fill="RoyalBlue")+xlim(0,3000) + labs(title="hg38 CGI region cg distance shift number = 50")
dev.off()



CGI <- fread("~/projects/CpGden/CpGIsland.hg38.txt")
setnames(CGI,c("chrom","start","end"))
CGI <- CGI[chrom==chr]
setkey(CGI,start)
CGI[,non.end := start-1]
CGI[,non.start := end+1]
non.CGI <- data.table(start=c(1,CGI$non.start),end=c(CGI$non.end,length(chr.seq)))
non.CGI[,id := c(1:nrow(non.CGI))]
non.CGI <- non.CGI[,.(pos=start:end),by= id]
non.CGI.distance <- merge(non.CGI,CG.info,by="pos")

pdf(file = "finalPlot/distancePlot/hg38-non-CGI-region-cg-distance.pdf")

ggplot(non.CGI.distance,aes(dis2))+geom_histogram(binwidth=1)+xlim(0,100) + labs(title="hg38 non CGI region cg distance shift number = 2")
ggplot(non.CGI.distance)+geom_density(aes(x=dis2),fill="RoyalBlue")+xlim(0,2000) + labs(title="hg38 non CGI region cg distance shift number = 2")

ggplot(non.CGI.distance,aes(dis5))+geom_histogram(binwidth=1)+xlim(0,2000) + labs(title="hg38 non CGI region cg distance shift number = 5")
ggplot(non.CGI.distance)+geom_density(aes(x=dis5),fill="RoyalBlue")+xlim(0,2000) + labs(title="hg38 non CGI region cg distance shift number = 5")
ggplot(non.CGI.distance,aes(dis10))+geom_histogram(binwidth=1)+xlim(0,5000) + labs(title="hg38 non CGI region cg distance shift number = 10")
ggplot(non.CGI.distance)+geom_density(aes(x=dis10),fill="RoyalBlue")+xlim(0,5000) + labs(title="hg38 non CGI region cg distance shift number = 10")
ggplot(non.CGI.distance,aes(dis15))+geom_histogram(binwidth=1)+xlim(0,5000) + labs(title="hg38 non CGI region cg distance shift number = 15")
ggplot(non.CGI.distance)+geom_density(aes(x=dis15),fill="RoyalBlue")+xlim(0,5000) + labs(title="hg38 non CGI region cg distance shift number = 15")
ggplot(non.CGI.distance,aes(dis20))+geom_histogram(binwidth=1)+xlim(0,5000) + labs(title="hg38 non CGI region cg distance shift number = 20")
ggplot(non.CGI.distance)+geom_density(aes(x=dis20),fill="RoyalBlue")+xlim(0,5000) + labs(title="hg38 non CGI region cg distance shift number = 20")
ggplot(non.CGI.distance,aes(dis25))+geom_histogram(binwidth=1)+xlim(0,5000) + labs(title="hg38 non CGI region cg distance shift number = 25") 
ggplot(non.CGI.distance)+geom_density(aes(x=dis25),fill="RoyalBlue")+xlim(0,5000) + labs(title="hg38 non CGI region cg distance shift number = 25")
ggplot(non.CGI.distance,aes(dis30))+geom_histogram(binwidth=1)+xlim(0,8000) + labs(title="hg38 non CGI region cg distance shift number = 30")
ggplot(non.CGI.distance)+geom_density(aes(x=dis30),fill="RoyalBlue")+xlim(0,8000) + labs(title="hg38 non CGI region cg distance shift number = 30")
ggplot(non.CGI.distance,aes(dis35))+geom_histogram(binwidth=1)+xlim(0,8000) + labs(title="hg38 non CGI region cg distance shift number = 35")
ggplot(non.CGI.distance)+geom_density(aes(x=dis35),fill="RoyalBlue")+xlim(0,8000) + labs(title="hg38 non CGI region cg distance shift number = 35")
ggplot(non.CGI.distance,aes(dis40))+geom_histogram(binwidth=1)+xlim(0,8000) + labs(title="hg38 non CGI region cg distance shift number = 40")
ggplot(non.CGI.distance)+geom_density(aes(x=dis40),fill="RoyalBlue")+xlim(0,8000) + labs(title="hg38 non CGI region cg distance shift number = 40")
ggplot(non.CGI.distance,aes(dis45))+geom_histogram(binwidth=1)+xlim(0,10000) + labs(title="hg38 non CGI region cg distance shift number = 45")
ggplot(non.CGI.distance)+geom_density(aes(x=dis45),fill="RoyalBlue")+xlim(0,10000) + labs(title="hg38 non CGI region cg distance shift number = 45")
ggplot(non.CGI.distance,aes(dis50))+geom_histogram(binwidth=1)+xlim(0,10000) + labs(title="hg38 non CGI region cg distance shift number = 50")
ggplot(non.CGI.distance)+geom_density(aes(x=dis50),fill="RoyalBlue")+xlim(0,10000) + labs(title="hg38 non CGI region cg distance shift number = 50")
dev.off()

pdf("finalPlot/distancePlot/hg38-chr1-noCGI-CGI-cg-distacnce-in-short-regrion-0-100-shitnu-2.pdf")
ggplot(CGI.distance,aes(dis2))+geom_histogram(binwidth=1)+xlim(0,100) + labs(title="hg38 CGI region cg distance shift number = 2")
ggplot(non.CGI.distance,aes(dis2))+geom_histogram(binwidth=1)+xlim(0,100) + labs(title="hg38 non CGI region cg distance shift number = 2")
ggplot(CGI.distance)+geom_density(aes(x=dis2),fill="RoyalBlue")+xlim(0,100) + labs(title="hg38 non CGI region cg distance shift number = 2")
ggplot(non.CGI.distance)+geom_density(aes(x=dis2),fill="RoyalBlue")+xlim(0,100) + labs(title="hg38 non CGI region cg distance shift number = 2")
dev.off()


#############################################
## pattern distance plot















## tg7
library(BSgenome.Tgondii.ToxoDB.7.0)
tg7 <- Tgondii
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(tg7) ){
  chr.seq <- DNAString(c(chr.seq,tg7[[tg.seq]]))
}
# wigPlot(chr.seq,pattern = "CG",wig.file = "finalPlot/tg7-chr1-CG-win2500-ps.wig",win = 1250,chrom = "chr1")
# wigPlot(chr.seq,pattern = "AT",wig.file = "finalPlot/tg7-chr1-AT-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "GC",wig.file = "finalPlot/tg7-chr1-GC-win2500-ps.wig",win = 1250,chrom = "chr1")
wigPlot(chr.seq,pattern = "TA",wig.file = "finalPlot/tg7-chr1-TA-win2500-ps.wig",win = 1250,chrom = "chr1")




## boxplot
##hg38
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38
hilbert.level <- 9
chr.seq <- hg38[[1]]
hg38.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
hg38.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
hg38.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
hg38.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
hg38.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")

##mm10
library(BSgenome.Mmusculus.UCSC.mm10)
mm10 <- BSgenome.Mmusculus.UCSC.mm10
hilbert.level <- 9
chr.seq <- mm10[[1]]
mm10.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
mm10.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
mm10.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
mm10.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
mm10.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")


## ce10
library(BSgenome.Celegans.UCSC.ce10)
ce10 <- BSgenome.Celegans.UCSC.ce10
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(ce10) ){
  chr.seq <- DNAString(c(chr.seq,ce10[[tg.seq]]))
}
ce10.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
ce10.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
ce10.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
ce10.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
ce10.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")


## danRer7
library(BSgenome.Drerio.UCSC.danRer7)
danRer7 <- BSgenome.Drerio.UCSC.danRer7
chr.seq <- danRer7[[1]]
hilbert.level <- 8
danRer7.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
danRer7.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
danRer7.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
danRer7.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
danRer7.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")


## dm6
library(BSgenome.Dmelanogaster.UCSC.dm6)
dm6 <- BSgenome.Dmelanogaster.UCSC.dm6
hilbert.level <- 9
chr.seq <- DNAString()
for(tg.seq in names(dm6)[!grepl("_", names(dm6))] ){
  chr.seq <- DNAString(c(chr.seq,dm6[[tg.seq]]))
}
dm6.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
dm6.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
dm6.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
dm6.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
dm6.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")

## ecoli
library(BSgenome.Ecoli.NCBI.20080805)
ecoli <- BSgenome.Ecoli.NCBI.20080805
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(ecoli) ){
  chr.seq <- DNAString(c(chr.seq,ecoli[[tg.seq]]))
}
ecoli.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
ecoli.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
ecoli.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
ecoli.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
ecoli.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")


## parkeri
library("Biostrings")
parkeri.path <- "~/projects/CpGden/parkeri.fasta"
parkeri <- readDNAStringSet(parkeri.path, format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=F)
hilbert.level <- 9
pkr <- DNAString()
for(tg.seq in c(1:100) ){
  # pkr <- DNAString(c(pkr,parkeri[tg.seq]))
  pkr <- DNAStringSet(paste(as.character(pkr),as.character(parkeri[tg.seq]),sep=""))
}
chr.seq <- pkr[[1]]
parkeri.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
parkeri.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
parkeri.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
parkeri.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
parkeri.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")


## sacCer3
library(BSgenome.Scerevisiae.UCSC.sacCer3)
saccer3<- Scerevisiae
hilbert.level <- 7
chr.seq <- DNAString()
for(tg.seq in names(saccer3) ){
  chr.seq <- DNAString(c(chr.seq,saccer3[[tg.seq]]))
}
sacCer3.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
sacCer3.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
sacCer3.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
sacCer3.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
sacCer3.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")

## taeGut1
library(BSgenome.Tguttata.UCSC.taeGut1)
taeGut1 <- BSgenome.Tguttata.UCSC.taeGut1
hilbert.level <- 8
chr.seq <- taeGut1[[1]]
taeGut1.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
taeGut1.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
taeGut1.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
taeGut1.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
taeGut1.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")

## tg7
library(BSgenome.Tgondii.ToxoDB.7.0)
tg7 <- Tgondii
hilbert.level <- 8
chr.seq <- DNAString()
for(tg.seq in names(tg7) ){
  chr.seq <- DNAString(c(chr.seq,tg7[[tg.seq]]))
}
tg7.CG.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "CG")
tg7.AT.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AT")
tg7.GC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "GC")
tg7.TA.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "TA")
tg7.AC.count <- binCount(chr.seq = chr.seq ,level = hilbert.level,pattern = "AC")
