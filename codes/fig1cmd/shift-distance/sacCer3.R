library(data.table)
library(HilbertCurve)
library(circlize)
library(ComplexHeatmap)
library(methyutils)
library(ggplot2)

## sacCer3
library(BSgenome.Scerevisiae.UCSC.sacCer3)
saccer3<- Scerevisiae
hilbert.level <- 7
chr.seq <- DNAString()
for(tg.seq in names(saccer3) ){
  chr.seq <- DNAString(c(chr.seq,saccer3[[tg.seq]]))
}

specie = "sacCer3"

##sacCer3 CG 
target.pattern="CG"
tp.pos <- matchPattern(target.pattern,chr.seq)
tp.pos <- start(tp.pos)
tp.shif <- data.frame(pos=tp.pos, 
                      sh2  = c(tp.pos[2:length(tp.pos)],tail(tp.pos,1)),
                      sh5  = c(tp.pos[5:length(tp.pos)],tail(tp.pos,4)),
                      sh10 = c(tp.pos[10:length(tp.pos)],tail(tp.pos,9)),
                      sh15 = c(tp.pos[15:length(tp.pos)],tail(tp.pos,14)),
                      sh20 = c(tp.pos[20:length(tp.pos)],tail(tp.pos,19)),
                      sh25 = c(tp.pos[25:length(tp.pos)],tail(tp.pos,24)),
                      sh30 = c(tp.pos[30:length(tp.pos)],tail(tp.pos,29)),
                      sh35 = c(tp.pos[35:length(tp.pos)],tail(tp.pos,34)),
                      sh40 = c(tp.pos[40:length(tp.pos)],tail(tp.pos,39)),
                      sh45 = c(tp.pos[45:length(tp.pos)],tail(tp.pos,44)),
                      sh50 = c(tp.pos[50:length(tp.pos)],tail(tp.pos,49)))


pdf(paste(specie,"-",target.pattern,"-distance-figure.pdf",sep=""))
theme_set(theme_bw(12))
ggplot(tp.shif,aes(sh20-pos))+geom_histogram(binwidth=1)+xlim(0,1300)+labs(title=paste(specie,target.pattern,"distance shift number = ",strsplit(names(tp.shif)[6],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh20-pos),fill="RoyalBlue")+xlim(0,1300) +labs(title=paste(specie,target.pattern,"distance shift number = ",strsplit(names(tp.shif)[6],"h")[[1]][2]))
dev.off()


pdf(paste(specie,"-",target.pattern,"-distance.pdf",sep=""))
ggplot(tp.shif,aes(sh2-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[2],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh2-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[2],"h")[[1]][2]))
ggplot(tp.shif,aes(sh5-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[3],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh5-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[3],"h")[[1]][2]))
ggplot(tp.shif,aes(sh10-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[4],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh10-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[4],"h")[[1]][2]))
ggplot(tp.shif,aes(sh15-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[5],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh15-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[5],"h")[[1]][2]))
ggplot(tp.shif,aes(sh20-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[6],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh20-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[6],"h")[[1]][2]))
ggplot(tp.shif,aes(sh25-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[7],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh25-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[7],"h")[[1]][2]))
ggplot(tp.shif,aes(sh30-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[8],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh30-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[8],"h")[[1]][2]))
ggplot(tp.shif,aes(sh35-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[9],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh35-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[9],"h")[[1]][2]))
dev.off()

##sacCer3 AT 
target.pattern="AT"
tp.pos <- matchPattern(target.pattern,chr.seq)
tp.pos <- start(tp.pos)
tp.shif <- data.frame(pos=tp.pos, 
                      sh2  = c(tp.pos[2:length(tp.pos)],tail(tp.pos,1)),
                      sh5  = c(tp.pos[5:length(tp.pos)],tail(tp.pos,4)),
                      sh10 = c(tp.pos[10:length(tp.pos)],tail(tp.pos,9)),
                      sh15 = c(tp.pos[15:length(tp.pos)],tail(tp.pos,14)),
                      sh20 = c(tp.pos[20:length(tp.pos)],tail(tp.pos,19)),
                      sh25 = c(tp.pos[25:length(tp.pos)],tail(tp.pos,24)),
                      sh30 = c(tp.pos[30:length(tp.pos)],tail(tp.pos,29)),
                      sh35 = c(tp.pos[35:length(tp.pos)],tail(tp.pos,34)),
                      sh40 = c(tp.pos[40:length(tp.pos)],tail(tp.pos,39)),
                      sh45 = c(tp.pos[45:length(tp.pos)],tail(tp.pos,44)),
                      sh50 = c(tp.pos[50:length(tp.pos)],tail(tp.pos,49)))

pdf(paste(specie,"-",target.pattern,"-distance-figure.pdf",sep=""))
theme_set(theme_bw(12))
ggplot(tp.shif,aes(sh20-pos))+geom_histogram(binwidth=1)+xlim(0,450)+labs(title=paste(specie,target.pattern,"distance shift number = ",strsplit(names(tp.shif)[6],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh20-pos),fill="RoyalBlue")+xlim(0,450) +labs(title=paste(specie,target.pattern,"distance shift number = ",strsplit(names(tp.shif)[6],"h")[[1]][2]))
dev.off()

pdf(paste(specie,"-",target.pattern,"-distance.pdf",sep=""))
ggplot(tp.shif,aes(sh2-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[2],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh2-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[2],"h")[[1]][2]))
ggplot(tp.shif,aes(sh5-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[3],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh5-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[3],"h")[[1]][2]))
ggplot(tp.shif,aes(sh10-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[4],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh10-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[4],"h")[[1]][2]))
ggplot(tp.shif,aes(sh15-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[5],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh15-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[5],"h")[[1]][2]))
ggplot(tp.shif,aes(sh20-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[6],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh20-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[6],"h")[[1]][2]))
ggplot(tp.shif,aes(sh25-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[7],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh25-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[7],"h")[[1]][2]))
ggplot(tp.shif,aes(sh30-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[8],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh30-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[8],"h")[[1]][2]))
ggplot(tp.shif,aes(sh35-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[9],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh35-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[9],"h")[[1]][2]))
dev.off()
##sacCer3 GC 
target.pattern="GC"
tp.pos <- matchPattern(target.pattern,chr.seq)
tp.pos <- start(tp.pos)
tp.shif <- data.frame(pos=tp.pos, 
                      sh2  = c(tp.pos[2:length(tp.pos)],tail(tp.pos,1)),
                      sh5  = c(tp.pos[5:length(tp.pos)],tail(tp.pos,4)),
                      sh10 = c(tp.pos[10:length(tp.pos)],tail(tp.pos,9)),
                      sh15 = c(tp.pos[15:length(tp.pos)],tail(tp.pos,14)),
                      sh20 = c(tp.pos[20:length(tp.pos)],tail(tp.pos,19)),
                      sh25 = c(tp.pos[25:length(tp.pos)],tail(tp.pos,24)),
                      sh30 = c(tp.pos[30:length(tp.pos)],tail(tp.pos,29)),
                      sh35 = c(tp.pos[35:length(tp.pos)],tail(tp.pos,34)),
                      sh40 = c(tp.pos[40:length(tp.pos)],tail(tp.pos,39)),
                      sh45 = c(tp.pos[45:length(tp.pos)],tail(tp.pos,44)),
                      sh50 = c(tp.pos[50:length(tp.pos)],tail(tp.pos,49)))

pdf(paste(specie,"-",target.pattern,"-distance.pdf",sep=""))
ggplot(tp.shif,aes(sh2-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[2],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh2-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[2],"h")[[1]][2]))
ggplot(tp.shif,aes(sh5-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[3],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh5-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[3],"h")[[1]][2]))
ggplot(tp.shif,aes(sh10-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[4],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh10-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[4],"h")[[1]][2]))
ggplot(tp.shif,aes(sh15-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[5],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh15-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[5],"h")[[1]][2]))
ggplot(tp.shif,aes(sh20-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[6],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh20-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[6],"h")[[1]][2]))
ggplot(tp.shif,aes(sh25-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[7],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh25-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[7],"h")[[1]][2]))
ggplot(tp.shif,aes(sh30-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[8],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh30-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[8],"h")[[1]][2]))
ggplot(tp.shif,aes(sh35-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[9],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh35-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[9],"h")[[1]][2]))
dev.off()
##sacCer3 TA 
target.pattern="TA"
tp.pos <- matchPattern(target.pattern,chr.seq)
tp.pos <- start(tp.pos)
tp.shif <- data.frame(pos=tp.pos, 
                      sh2  = c(tp.pos[2:length(tp.pos)],tail(tp.pos,1)),
                      sh5  = c(tp.pos[5:length(tp.pos)],tail(tp.pos,4)),
                      sh10 = c(tp.pos[10:length(tp.pos)],tail(tp.pos,9)),
                      sh15 = c(tp.pos[15:length(tp.pos)],tail(tp.pos,14)),
                      sh20 = c(tp.pos[20:length(tp.pos)],tail(tp.pos,19)),
                      sh25 = c(tp.pos[25:length(tp.pos)],tail(tp.pos,24)),
                      sh30 = c(tp.pos[30:length(tp.pos)],tail(tp.pos,29)),
                      sh35 = c(tp.pos[35:length(tp.pos)],tail(tp.pos,34)),
                      sh40 = c(tp.pos[40:length(tp.pos)],tail(tp.pos,39)),
                      sh45 = c(tp.pos[45:length(tp.pos)],tail(tp.pos,44)),
                      sh50 = c(tp.pos[50:length(tp.pos)],tail(tp.pos,49)))

pdf(paste(specie,"-",target.pattern,"-distance.pdf",sep=""))
ggplot(tp.shif,aes(sh2-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[2],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh2-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[2],"h")[[1]][2]))
ggplot(tp.shif,aes(sh5-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[3],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh5-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[3],"h")[[1]][2]))
ggplot(tp.shif,aes(sh10-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[4],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh10-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[4],"h")[[1]][2]))
ggplot(tp.shif,aes(sh15-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[5],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh15-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[5],"h")[[1]][2]))
ggplot(tp.shif,aes(sh20-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[6],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh20-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[6],"h")[[1]][2]))
ggplot(tp.shif,aes(sh25-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[7],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh25-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[7],"h")[[1]][2]))
ggplot(tp.shif,aes(sh30-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[8],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh30-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[8],"h")[[1]][2]))
ggplot(tp.shif,aes(sh35-pos))+geom_histogram(binwidth=1)+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999))+labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[9],"h")[[1]][2]))
ggplot(tp.shif)+geom_density(aes(sh35-pos),fill="RoyalBlue")+xlim(0,quantile(with(tp.shif,sh5-pos),0.9999)) +labs(title=paste(specie,"cg distance shift number = ",strsplit(names(tp.shif)[9],"h")[[1]][2]))
dev.off()

