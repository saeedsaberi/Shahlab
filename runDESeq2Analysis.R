library(DESeq2)
library("RColorBrewer")
library("gplots")
library("vsn")



args <- commandArgs(trailingOnly = T)

wdir=as.character(args[1])
setwd(wdir)
inputDir <- wdir
col <-as.integer(args[3])
outputDir <- as.character(args[2])

#col <- "V11"

#inputDir <- "/Users/ssaberi/Documents/DLBCL/counts"
#listDir <- "/Users/ssaberi/Documents/DLBCL/counts"
#outputDir <- "/Users/ssaberi/Documents/DLBCL/counts/DESeq2CountResults/"

categories <- read.table(paste0(wdir,"categories.txt"),sep="\t")
myfun <- function(x){
  return(gsub("_","",x))
}
categoryID <- lapply(categories$V1,myfun)
categories$V1 <- categoryID

myfun <- function(x){
  fls<- list.files(inputDir)
  if (is.null(x) || is.na(x) ){
    return(NA)
    print(x) 
  }
  tmp=strsplit(x,"/")[[1]]
  x <- tmp[length(tmp)]
  if (length(grep(x,fls))>0){
    x<- paste0(x,"*")
    dir<- paste0(inputDir,"/")
    x <- paste0(dir,x)
    ret=Sys.glob(x)

    return(ret[1]) 
    }else {
      return(NA)
      
    }
  }

countFile1Temp <- lapply(categories$V1[categories[,col]==levels(categories[,col])[1]],myfun)
countFile1Temp <- countFile1Temp[!is.na(countFile1Temp)]
countFile1Temp <- countFile1Temp[!is.null(countFile1Temp)]
countFile2Temp <- lapply(categories$V1[categories[,col]==levels(categories[,col])[2]],myfun)
countFile2Temp <- countFile2Temp[!is.na(countFile2Temp)]
countFile2Temp <- countFile2Temp[!is.null(countFile2Temp)]

countFiles <- c(t(countFile1Temp), t(countFile2Temp))
sampleFiles <- c()
for(i in 1:length(countFiles)){
  sampleFiles<-rbind(sampleFiles,countFiles[i])
}

sampleConditions<- c()
for(i in 1:length(countFile1Temp)){
  sampleConditions<-rbind(sampleConditions,levels(categories[,col])[1])
}
for(i in 1:length(countFile2Temp)){
  sampleConditions<-rbind(sampleConditions,levels(categories[,col])[2])
}
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, sampleConditions=sampleConditions)
#View(sampleTable)


ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory =inputDir, design = ~ sampleConditions)
#print(ddsHTSeq)

dds <- DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)

pdf(paste0(outputDir,"DESeq2_MA.pdf"))
plotMA(dds,ylim=c(-2,2),main="DESeq2.pdf")
dev.off()
resOrdered <- res[order(res$padj), ]

tmp<-subset(resOrdered,resOrdered$padj<.01)
upgenes<-subset(tmp,tmp$log2FoldChange>= .4)
downgenes<-subset(tmp,tmp$log2FoldChange<= -.4)

write.table(resOrdered,paste0(outputDir,"simple-res.txt"),sep="\t")
write.table(upgenes,paste0(outputDir,"simple-res-upgenes.txt"),sep="\t")
write.table(downgenes,paste0(outputDir,"simple-res-downgenes.txt"),sep="\t")


vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

pdf(paste0(outputDir,"DESeq2_vsd_variance.pdf"))
par(mfrow=c(1,2))
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1), ylim = c(0,2.5))
meanSdPlot(assay(vsd[notAllZero,]), ylim = c(0,2.5))

title(main="shifted(dds) log_2")

meanSdPlot(assay(vsd[notAllZero,]), ylim = c(0,2.5))
title(main="vsd variance")
dev.off()

pdf(paste0(outputDir,"DESeq2_heatmap1.pdf"))
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "RdYlBu"))(50)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))

dev.off()

pdf(paste0(outputDir,"DESeq2_heatmap3.pdf"))
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))

dev.off()



#distsv <- dist(t(assay(vsd)))
distsv <-as.dist(1-cor(assay(vsd), method="spearman"))
mat <- as.matrix(distsv)
hc <- hclust(distsv, method="complete")
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(sampleConditions,sampleFiles , sep=" : "))


pdf(paste0(outputDir,"deseq2_heatmaps_samplebysample.pdf"))
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, 
          col = rev(hmcol), margin=c(13, 13),
          distfun=function(c) as.dist(1 - c), trace="none",
          RowSideColors=sampleConditions)
legend("topright",      
               legend = unique(sampleConditions),
               col = unique(as.numeric(sampleConditions)), 
               lty= 1,             
               lwd = 5,           
               cex=.7 )

dev.off()


pdf(paste0(outputDir,"deseq2_pca.pdf"))
print(plotPCA(vsd, intgroup=c("sampleConditions")))
dev.off()



### new
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)

tab <- table(initial = results(dds)$padj < .1,
             cleaned = results(ddsClean)$padj < .1)
addmargins(tab)
write.csv(as.data.frame(tab),file=paste0(outputDir,
                       "sim_condition_treated_results_cleaned_summary_deseq2.csv"))
resClean <- results(ddsClean)
write.csv(as.data.frame(resClean),file=paste0(outputDir,
                                      "sim_condition_treated_results_cleaned_deseq2.csv"))

#filtering threashold
pdf(paste0(outputDir,"deseq2_filtering_treshold.pdf"))
attr(res,"filterThreshold")
#     10%
#91.48005
plot(attr(res,"filterNumRej"),type="b", ylab="number of rejections")
dev.off()

W <- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
pdf(paste0(outputDir,"deseq2_cooksdist.pdf"))
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic",
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))
dev.off()

pdf(paste0(outputDir,"deseq2_indep_filt.pdf"))
plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
dev.off()

pdf(paste0(outputDir,"barplots-cleaned.pdf"))
use <- res$baseMean > attr(res,"filterThreshold")
table(use)
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c('do not pass'="khaki", 'pass'="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))
dev.off()


resFilt <- res[use & !is.na(res$pvalue),]
orderInPlot <- order(resFilt$pvalue)
showInPlot <- (resFilt$pvalue[orderInPlot] <= 0.08)
alpha <- 0.1

pdf(paste0(outputDir,"barplots-cleaned.pdf"))
plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],
     pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))
abline(a=0, b=alpha/length(resFilt$pvalue), col="red3", lwd=2)
dev.off()