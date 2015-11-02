

#load edgeR library
#load biomaRt library
#load the deseq library
library(DESeq)
library("gplots")
library("RColorBrewer")
library("vsn")

countsTable <- read.table("~/Documents/DLBCL/DLBCL-ALL.txt",sep='\t')
categories <- read.table("~/Documents/DLBCL/categories.txt",sep="\t")
myfun <- function(x){
return(gsub("_","",x))
}
categoryID <- lapply(categories$V1,myfun)
categories$V1 <- categoryID



vec <- vector(mode="integer", length=length(countsTable))
for (i in 1:length(countsTable)){
  vec[i]=categories$V11[which(categories$V1==colnames(countsTable)[i])]
}


counttry<- countsTable[vec<3]
vec=vec[vec<3]

counttry <- (counttry[,colnames(counttry)!="DLC0260"])
counttry <- (counttry[,colnames(counttry)!="DLC0262"])
counttry <- (counttry[,colnames(counttry)!="DLC0169"])
counttry <- (counttry[,colnames(counttry)!="DLC0224"])
counttry <- (counttry[,colnames(counttry)!="DLC0243"])
counttry <- (counttry[,colnames(counttry)!="DLC0219"])

nsamp= length(counttry)
vec <- vector(mode="integer", length=length(counttry[1:nsamp]))
for (i in 1:length(counttry[1:nsamp])){
  vec[i]=categories$V11[which(categories$V1==colnames(counttry)[i])]
}

conds <- factor( c(vec))
cds <- newCountDataSet( counttry[,1:nsamp], conds )
cds <- DESeq::estimateSizeFactors( cds )
sizeFactors( cds )
cds <- DESeq::estimateDispersions( cds )

pdf('~/Documents/DLBCL/DR.pdf')
DESeq::plotDispEsts( cds )
dev.off()



pdf("~/Documents/DLBCL/plotMA.pdf")
res = DESeq::nbinomTest( cds, 1,2 )
DESeq::plotMA(res)
dev.off()


ncu = counts( cds, normalized=TRUE )[ , conditions(cds)==1 ]

jpeg("MA_ABC.jpg")
DESeq::plotMA(data.frame(baseMean = rowMeans(ncu),log2FoldChange = log2( ncu[,2] / ncu[,1] )),
        col = "black")
dev.off()


hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")

jpeg("hist_pval.jpg")
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()

resSig = res[ res$padj < 0.1, ]
head(resSig)
write.csv( res, file="72hour_unstgVstg_all.csv" )
write.csv( resSig, file="72hour_unstgVstg_only_0.01FDR.csv" )


rs = rowSums ( counts ( cds))

theta = 0.4

use = (rs > quantile(rs, probs=theta))

table(use)
cdsFilt = cds[ use, ]
resFilt = DESeq::nbinomTest( cdsFilt,1, 2 )
tabFull
tabFilt

jpeg("pval_bot40.jpg")

plot(rank(rs)/length(rs), -log10(res$pval), pch=16, cex=0.45)

dev.off()

h1 = hist(res$pval,breaks=50,plot=FALSE)

h2 = hist(resFilt$pval,breaks=50,plot=FALSE)

colori <- c('do not pass'="khaki", 'pass'='powderblue')

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori,
          + space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))
dev.copy(png,’hist_filthist.png’)
dev.off()



cdsBlind = estimateDispersions( cds, method="blind" )
vsd = varianceStabilizingTransformation( cdsBlind )

par(mfrow=c(1,2))
notAllZero = (rowSums(counts(cds))>0)
meanSdPlot(log2(counts(cds)[notAllZero, ] + 1), ylim = c(0,2.5))
meanSdPlot(vsd[notAllZero, ], ylim = c(0,2.5))

mod_lfc = (rowMeans( exprs(vsd)[, conditions(cds)==2, drop=FALSE] ) – rowMeans( exprs(vsd)[, conditions(cds)==1, drop=FALSE] ))
lfc = res$log2FoldChange
table(lfc[!is.finite(lfc)], useNA="always")

logdecade = 1 + round( log10( 1+rowMeans(counts(cdsBlind, normalized=TRUE)) ) )
lfccol = colorRampPalette( c( "gray", "blue" ) )(6)[logdecade]
ymax = 4.5
plot( pmax(-ymax, pmin(ymax, lfc)), mod_lfc,
        xlab = "ordinary log-ratio", ylab = "moderated log-ratio",
        cex=0.45, asp=1, col = lfccol,
        pch = ifelse(lfc<(-ymax), 60, ifelse(lfc>ymax, 62, 16)))

abline( a=0, b=1, col="red3")


cdsBlind = estimateDispersions(cds,method="blind")
vsd = varianceStabilizingTransformation(cdsBlind)
select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsd)[select,], col = hmcol, trace="none", margin=c(10, 6))
