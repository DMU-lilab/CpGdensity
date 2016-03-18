library("data.table")

overlap.region <- read.table("~/cgDensity/IcmpRegion/Sw2500Region/data2/bonemarrow/regionValue/overlap/region_mtbr_value.txt",header=T)

overlap.region <- overlap.region[overlap.region$bonemarrow_cgCover > 20,]

library("matrixStats") 

s <- grep("score",colnames(overlap.region ))
ts.score <- overlap.region[,s]
ts.score.matrix <- as.matrix(ts.score)
ts.score.sd <- rowSds(ts.score.matrix)
overlap.region$sd  <- ts.score.sd
overlap.region <- overlap.region[!is.na(overlap.region$sd),]
overlap.region.sdorder <- overlap.region[order(overlap.region$sd, decreasing = T),]



length(overlap.region.sdorder$sd )  #9797
sum(overlap.region.sdorder$sd >= 0.2) # 91 0.93%


##############
gap.value <- read.table("~/cgDensity/IcmpRegion/Sw2500Region/data3/bonemarrow/regionValue/gap/region_mtbr_value.txt",header = T)
gap.value <- gap.value[gap.value$bonemarrow_cgCover > 5,]
#gap.value.bone <- gap.value[,c(1:9)]
gap.value.bone <- gap.value
#dt.gap.value.bone <- as.data.table(gap.value.bone)
#bone.cg.freq <- as.data.table(table(dt.gap.value.bone$bonemarrow_cgCover))
gap.value.bone$length <- gap.value.bone$end - gap.value.bone$start
gap.value.bone$cgDen <- gap.value.bone$bonemarrow_cgCover / gap.value.bone$length
gap.value.bone <- gap.value.bone[gap.value.bone$cgDen > 0.003,]

s <- grep("score",colnames(gap.value.bone ))
ts.score <- gap.value.bone[,s]
ts.score.matrix <- as.matrix(ts.score)
ts.score.sd <- rowSds(ts.score.matrix)
gap.value.bone$sd  <- ts.score.sd
gap.value.bone <- gap.value.bone[!is.na(gap.value.bone$sd),]
gap.value.bone.sdorder <- gap.value.bone[order(gap.value.bone$sd, decreasing = T),]


################
tissues <- c("bonemarrow","cerebellum","skin","stomach","colon","cortex","heart","intestine","kidney","liver","lung","olfactorybulb","pancreas","placenta","spleen","thymus","uterus")
chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
#############################################
a <- data.table(chrom = "chr1",start = c(2,10,20),end = c(5,15,25),id = c(1,2,3)) 
b <- a[,.(chrom = chrom , pos = c(start:end)),by = id] 
value <- data.table(pos = c(1:30),value =c(1:30)) 
c <-  merge(b,value,by = "pos")
c[,median(value),by = id]

#################overlap
overlap.region.bone <- overlap.region[,c(1:3)]
chrs.icmp.list <- list()
for (chr in chrs){
  overlap.region.bone.chr <- overlap.region.bone[overlap.region.bone$chrom == chr,]
  overlap.region.bone.chr$id <- 1:nrow(overlap.region.bone.chr)
  dt.overlap.region.bone.chr <- as.data.table(overlap.region.bone.chr)
  dt.overlap.region.bone.chr.split <- dt.overlap.region.bone.chr[,.(chrom = chrom,pos = c(start : end)), by = id]
  
  tissues.icmp.list <- list()
  for (ts in tissues){
    message(chr,"  ", ts, " is running ",date())
    load(paste("~/cgDensity/IcmpRegion/Sw2500Region/ICMP/", ts, "/", chr,".Rdata",sep = ""))
    ts.value <- data.table(pos = 1:length(icmp), icmp = icmp)
    dt.overlap.region.ts.value <- merge(dt.overlap.region.bone.chr.split, ts.value, by = "pos")
    ts.icmp <- dt.overlap.region.ts.value[, median(icmp), by = id]
    tissues.icmp.list[[ts]] <- ts.icmp$V1

  }
  tissues.icmp <- do.call(cbind,tissues.icmp.list)
  dt.overlap.region.chr.icmp <- cbind(dt.overlap.region.bone.chr,tissues.icmp)
  chrs.icmp.list[[chr]] <- dt.overlap.region.chr.icmp
}
overlap.region.bone.icmp <- do.call(rbind, chrs.icmp.list)
df.overlap.region.bone.icmp <- as.data.frame(overlap.region.bone.icmp)


test <- df.overlap.region.bone.icmp[,c(5:21)]
test <- test[,c(1:13,15:17)]
overlap.func <- function(x){
  c <- c()
  for (a in x){
    b <- a < 0 & a > -0.3
    c <- c(c,b)
    s <- sum(c)
  }
  return(s)  
}
test4 <- apply(test,1,overlap.func)

df.overlap.region.bone.icmp$icmpmaker <- test3
df.overlap.region.bone.icmp.un <- df.overlap.region.bone.icmp[df.overlap.region.bone.icmp$icmpmaker == 0,] ### sum(no.region)383 9.7% sum(gapValueRepeat) 3529 90.3%
df.overlap.region.bone.icmp.on <- df.overlap.region.bone.icmp[df.overlap.region.bone.icmp$icmpmaker > 0,] ### sum(no.region) 811 13.8% sum(gapValueRepeat) 5061 86.2%


###############gap region
test3 <- apply(test,1,overlap.func)
gap.value.bone$icmpmaker <- test3

#############
### pie
overlap.icmpmaker <- as.data.frame(table(df.overlap.region.bone.icmp$icmpmaker))
overlap.icmpmaker$maker <- paste("(",overlap.icmpmaker$Var1,")",round(overlap.icmpmaker$Freq / sum(overlap.icmpmaker$Freq)* 100,1), "%",sep="") 
overlap.icmpmaker$summaker <- c(0,rep(1:3,each = 4),4,4)
dt.overlap.icmpmaker <- data.table(Freq = overlap.icmpmaker$Freq, summaker = overlap.icmpmaker$summaker) 
dt.overlap.icmpmaker.new <- dt.overlap.icmpmaker[,.(sum = sum(Freq)),by = summaker]
dt.overlap.icmpmaker.new$Type <- c("0","1-4","5-8","9-12","13-14") 

colours = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
library(ggplot2)

percent_str <- paste(round(dt.overlap.icmpmaker.new$sum/sum(dt.overlap.icmpmaker.new$sum) * 100,1), "%", sep="")
dt.overlap.icmpmaker.new$Percentage <- percent_str
pie <-
  ggplot(dt.overlap.icmpmaker.new, aes(x = "" ,y = sum, fill = Type)) +  geom_bar(width = 3, stat = "identity")  + coord_polar("y") + xlab('') + ylab('') + labs(fill="Types") + scale_fill_manual(values = colours) + scale_y_continuous(breaks = NULL) + publication.theme 

publication.theme <- theme_bw(15) + theme(axis.title.y=element_text(vjust=1.7), axis.title.x=element_text(vjust=-0.1), text= element_text(size = 24, face = "bold"), axis.line = element_line(colour = "black", size = 0), panel.border = element_rect(colour = "black", size = 0, fill = NA), panel.background = element_rect(fill = "white", size = 5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))
                                   
library("RColorBrewer")
display.brewer.pal
brewer.pal(9,"Set1")