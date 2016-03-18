### 2015/12/6 write the code for selecting gap and overlap region
library("data.table")
library("ggplot2")
library("matrixStats") 
gap.value <- read.table("~/cgDensity/IcmpRegion/Sw2500Region/data/bonemarrow/regionValue/gap/region_mtbr_value.txt",header=T)
gap.value <- gap.value[gap.value$bonemarrow_cgCover > 5,]
#gap.value.bone <- gap.value[,c(1:9)]
gap.value.bone <- gap.value
#dt.gap.value.bone <- as.data.table(gap.value.bone)
#bone.cg.freq <- as.data.table(table(dt.gap.value.bone$bonemarrow_cgCover))
gap.value.bone$length <- gap.value.bone$end - gap.value.bone$start
gap.value.bone$cgDen <- gap.value.bone$bonemarrow_cgCover / gap.value.bone$length
gap.value.bone <- gap.value.bone[gap.value.bone$cgDen > 0.003,]

pdf("~/cgDensity/IcmpRegion/Sw2500Region/figure/gap/bonemarrow/bonemarrow_cgCover_5_50.pdf")
ggplot(gap.value.bone,aes(x=bonemarrow_cgCover)) + geom_histogram(binwidth=1) + xlim(5,50)
dev.off()

pdf("~/cgDensity/IcmpRegion/Sw2500Region/figure/gap/bonemarrow/bonemarrow_regionLength.pdf")
ggplot(gap.value.bone,aes(x=length)) + geom_histogram(binwidth=100) + xlim(0,10000)
dev.off()

tissues <- c("spleen","thymus","bonemarrow","colon","intestine","pancreas","stomach","liver","heart","lung","kidney","uterus","skin","cerebellum","cortex","olfactorybulb","placenta")
for (ts in tissues){
  pdf(paste("~/cgDensity/IcmpRegion/Sw2500Region/figure/gap/bonemarrow/tsScoreHistogram/",ts,".pdf",sep=""))
  a <- paste(ts,"_score",sep="")
  ggplot(gap.value.bone,aes(x= a)) + geom_histogram(binwidth=0.02) + xlim(0,1)
  dev.off()
  }

library("matrixStats") 
s <- grep("score",colnames(orderrvd))
score <- orderrvd[,s,with=FALSE]  
scorematrix <- as.matrix(score) 
sd <- rowSds(scorematrix)
orderrvd$sd <- sd

s <- grep("score",colnames(gap.value.bone ))
ts.score <- gap.value.bone[,s]
ts.score.matrix <- as.matrix(ts.score)
ts.score.sd <- rowSds(ts.score.matrix)
gap.value.bone$sd  <- ts.score.sd
gap.value.bone <- gap.value.bone[!is.na(gap.value.bone$sd),]
gap.value.bone.sdorder <- gap.value.bone[order(gap.value.bone$sd, decreasing = T),]

write.table(gap.value.bone.sdorder,"~/cgDensity/IcmpRegion/Sw2500Region/data/bonemarrow/regionValue/gap/bonemarrow_select_gapRegionValue.txt",col.names=T,row.names=F,sep = "\t")


### find the interesting gene in the gap region 
library(IRanges)

geneValue <- read.table("/home/fsch/cgDensity/mm9genepromoter/intesetinggene_valueSd0.2.txt",header=T)
geneValue <- data.frame(g:weneNum=c(1:126),geneValue) 

chrom <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX")

chr.Region.list <- list()
for (chr in chrom){
	chr.gap.value.bone <- gap.value.bone[gap.value.bone$chrom == chr,]
	gapRegion <- chr.gap.value.bone[,1:3]
	
	chr.gene.value <- geneValue[geneValue$chrom == chr,]
	##findOverlaps
	query <- IRanges(chr.gene.value$regionStart, chr.gene.value$regionEnd)
	subject <- IRanges(gapRegion$start, gapRegion$end)
	overlap <- findOverlaps(query,subject)
	
	chr.gene.value.overlap <- chr.gene.value[overlap@queryHits,]
	gapRegion.overlap <- gapRegion[overlap@subjectHits,]
	gapInterest.Gene <- cbind(gapRegion.overlap, chr.gene.value.overlap)
	chr.Region.list[[chr]] <- gapInterest.Gene
}    
chr.Region <- do.call(rbind, chr.Region.list)






################# overlap region 
library("data.table")

overlap.region <- read.table("~/cgDensity/IcmpRegion/Sw2500Region/data/bonemarrow/regionValue/overlapValue/region_mtbr_value.txt",header=T)
mm9.repeat <- fread("~/ws03-backup/fsch/cgDensity/RepeatMasker/repeatmm9",header=T,sep="\t")

overlap.region <- overlap.region[overlap.region$bonemarrow_cgCover > 20,]
overlap.region.bone <- overlap.region[,c(1:9)]

chrom <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX")

no.region <- c()
overlapValueRepeatlist <- list()
for (chr in chrom) {

        overlap.region.bone.chr <- overlap.region.bone[overlap.region.bone$chrom==chr,]
        #load(paste("/home/fsch/ws03-backup/fsch/cgDensity/RepeatMasker/Rdata/",chr,".Rdata",sep=""))
        mm9.repeat.chr <- mm9.repeat[mm9.repeat$genoName == chr,]
	query <- IRanges(overlap.region.bone.chr$start,overlap.region.bone.chr$end)
        subject <- IRanges(mm9.repeat.chr$genoStart,mm9.repeat.chr$genoEnd)
        findoverlap <- findOverlaps(query,subject)
        findRepeat <- data.table(findoverlap@queryHits,findoverlap@subjectHits,mm9.repeat.chr$repClass[findoverlap@subjectHits])       
        #findRepeatNew <- findRepeat[,list(class=paste(unique(V3),collapse=",")),by=list(V1)]
	findRepeatNew <- findRepeat[,list(class=paste(unique(V3))),by=list(V1)]
	#overlap.region.bone.chr$RepeatClass <- "no"
        #overlap.region.bone.chr$RepeatClass[findRepeatNew$V1] <- findRepeatNew$class
	overlapValueRepeatlist[[chr]] <- findRepeatNew
	no.num <- length(overlap.region.bone.chr$chrom) - length(unique(findRepeat$V1))
	no.region <- c(no.region,no.num)
}
Repeat.value <- do.call(rbind,overlapValueRepeatlist)
df.repeat <- as.data.frame(table(Repeat.value$class))
colnames(df.repeat) <- c("class","num")
df.no <- data.frame(class = "NO", num = sum(no.region))

df.class.num <- rbind(df.repeat,df.no)
df.class.num.order <- df.class.num[order(df.class.num$num),]
Others.num <- sum(df.class.num.order$num[1:9])
df.class <- df.class.num.order[10:16,]
df.others <- data.frame(class = "Others", num = Others.num)
df.class.new <- rbind(df.class,df.others)
df.class.new$class.percent <- class.percent <- paste(df.class.new$class,round(df.class.new$num / sum(df.class.new$num) * 100,1), "%",sep="")  

library(MASS)

pie(df.class.new$num,labels = df.class.new$class.percent)



