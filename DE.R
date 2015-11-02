#limma differential expression analysis
#load Rsubread library
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite(c("GenomicFeatures", "AnnotationDbi","DESeq2","DESeq","edgeR","Rsamtools","Rsubread","biomaRt"))

library(Rsubread)
#load edgeR library
library(edgeR)
#load biomaRt library
library(biomaRt)
#load the deseq library
library(DESeq)
library(DESeq2)
library("Rsamtools")


files<-read.csv("bamfiles.txt",header=F)
files <- data.frame(files)
bamFls<-as.vector(files$V1)

count.mat <- featureCounts(bamFls,
                           annot.ext="/share/lustre/reference/genomes/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf",
                           isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id", useMetaFeatures=TRUE,isPairedEnd=TRUE)
write.table(count.mat,file="/home/ssaberi/DLBCL/DLBCL-322.bams.count",sep="\t")

#reading the samples
samples<- read.csv("samples.txt",header=F)
samples <- data.frame(samples)
samples<-as.vector(samples$V1)


count.col<-samples
colnames(count.mat$counts)<-count.col






#choose the homosapiens datamart
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#read in the coding genes data set from the nature paper
nature.counts<-read.delim(
       file="/share/lustre/kshumansky/Alignments/Huntsman/FF/wtss/diffExpr/EGA_cellLines/140625_Klijn_counts_coding.txt",
       header=T, stringsAsFactors=FALSE)




