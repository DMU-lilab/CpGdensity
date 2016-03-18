### 2015/12/10 write the code for gap region with gene 

gene.ts <- read.table("/home/fsch/ws03-backup/fsch/mm9gene_promoter/resultValue/promoter_2000/intesetinggene_tissues_value.txt",sep = "\t", header = T)
gene.num <- table(gene.ts$geneSymbols)
gene.num <- as.data.frame(gene.num)
colnames(gene.num) <- c("geneSymbols", "geneNum")

gene.ts <- merge(gene.ts, gene.num, by = "geneSymbols")

gene.ts.2 <- gene.ts[gene.ts$geneNum == 2,] #gene num 70 / 2 = 35
gene.ts.3 <- gene.ts[gene.ts$geneNum == 3,] # gene num 93 / 3 = 31
gene.ts.4 <- gene.ts[gene.ts$geneNum == 4,] # gene num 84 / 4 = 21
gene.ts.5 <- gene.ts[gene.ts$geneNum == 5,] # gene num 30 / 5 = 6
gene.ts.6 <- gene.ts[gene.ts$geneNum == 6,] # gene num 42 / 6 = 7
gene.ts.7 <- gene.ts[gene.ts$geneNum == 7,] # gene num 21 / 7 = 3

###

s <- grep("score",colnames(gene.ts.2))
gene.ts.score.2 <- gene.ts.2[,s]
gene.ts.score.2$geneSymbols <- gene.ts.2$geneSymbols
gene.ts.score.2$sd <- gene.ts.2$sd 
rownames(gene.ts.score.2) <- c(1:70)
gene.ts.score.2 <- gene.ts.score.2[-c(23,24,35,36,39,40,69,70),]

gene.ts.score.2.diff <- gene.ts.score.2[gene.ts.score.2$sd > 0.2,]
##################

#################
plot.data1 <- gene.ts.score.2.diff[,1:17]
plot.data1.scale <- scale(plot.data1)
plot.data1.scale.kmeans <- kmeans(plot.data1.scale,4)

gene.ts.score.2.diff.cluster <- cbind(gene.ts.score.2.diff,cluster = plot.data1.scale.kmeans$cluster)
##########################
gene.cluster <- gene.ts.score.2.diff.cluster[,c(18,20)]
gene.ts.cluster <- merge(gene.ts.score.2,gene.cluster,by = "geneSymbols")
gene.ts.cluster.order <- gene.ts.cluster[order(gene.ts.cluster$cluster,gene.ts.cluster$geneSymbols),]
mm9geneheatorder <- gene.ts.cluster.order[,c(15,17,2,4,7,12,16,9,6,10,8,18,14,3,5,11,13)] 
#################################################
gene.ts.score.2.diff.cluster.order <- gene.ts.score.2.diff.cluster[order(gene.ts.score.2.diff.cluster$cluster),]


mm9geneheatorder <- gene.ts.score.2.diff.cluster.order[,c(14,16,1,3,6,11,15,8,5,9,7,17,13,2,4,10,12)] 
colnames(mm9geneheatorder) <- c("spleen","thymus","bonemarrow","colon","intestine","pancreas","stomach","liver","heart","lung","kidney","uterus","skin","cerebellum","cortex","olfactorybulb","placenta")


mat <- data.matrix(mm9geneheatorder[,1:17])
row.names(mat) <- c(1:nrow(mat)) 
m <- melt(mat) 
#png("/home/cfs/project/cgDen/mm9gene/promoter_2000/intesteringGene.png")
g <- ggplot(m,aes(x=Var2,y=Var1,fill=value)) + xlab('tissues') + ylab('genes')


g + geom_tile() +scale_fill_gradient2(low="green",mid="black",high="red",midpoint=0.55,limits=c(0,1)) +scale_y_continuous(expand = c(0, 0),labels=gene.ts.cluster.order$geneSymbols,breaks=1:length(gene.ts.cluster.order$geneSymbols))+ theme(axis.text.y=element_text(size=5)) + theme(axis.text.x=element_text(size=9,angle=90)) + publication.theme

publication.theme <- theme_bw(15) + theme(axis.title.y=element_text(vjust=1.7), axis.title.x=element_text(vjust=-0.1), text= element_text(size = 24, face = "bold"), axis.line = element_line(colour = "black", size = 0.1), panel.border = element_rect(colour = "black", size = 1, fill = NA), panel.background = element_rect(fill = "white", size = 5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

mm9.gene.1 <- mm9geneheatorder[seq(1,61,by = 2),]
colnames(mm9.gene.1) <- 1:17
mm9.gene.2 <- mm9geneheatorder[seq(2,62,by = 2),]
colnames(mm9.gene.2) <- 1:17
mm9.gene.combine <- cbind(mm9.gene.1,mm9.gene.2)
mm9.gene.combine.order <-  mm9.gene.combine[,order(as.numeric(colnames((mm9.gene.combine))))]
mat <- data.matrix(mm9.gene.combine.order)
colnames(mat) <- c("spleen1","spleen2","thymus1","thymus2","bonemarrow1","bonemarrow2","colon1","colon2","intestine1","intestine2","pancreas1","pancreas2","stomach1","stomach2","liver1","liver2","heart1","heart2","lung1","lung2","kidney1","kidney2","uterus1","uterus2","skin1","skin2","cerebellum1","cerebellum2","cortex1","cortex2","olfactorybulb1","olfactorybulb2","placenta1","placenta2")
rownames(mat) <- c(1:nrow(mm9.gene.combine.order)) 
m <- melt(mat) 
#png("/home/cfs/project/cgDen/mm9gene/promoter_2000/intesteringGene.png")
g <- ggplot(m,aes(x=Var2,y=Var1,fill=value)) + xlab('tissues') + ylab('genes')


g + geom_tile() +scale_fill_gradient2(low="green",mid="black",high="red",midpoint=0.55,limits=c(0,1)) +scale_y_continuous(expand = c(0, 0),labels=unique(gene.ts.cluster.order$geneSymbols),breaks=1:length(unique(gene.ts.cluster.order$geneSymbols)))+ theme(axis.text.y=element_text(size=5)) + theme(axis.text.x=element_text(size=9,angle=90))
