
#load biomaRt library
library(biomaRt)
#load the deseq library
library(DESeq2)


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

counttry <- countsTable


counttry <- (counttry[,colnames(counttry)!="DLC0260"])
counttry <- (counttry[,colnames(counttry)!="DLC0262"])
counttry <- (counttry[,colnames(counttry)!="DLC0169"])
counttry <- (counttry[,colnames(counttry)!="DLC0224"])
counttry <- (counttry[,colnames(counttry)!="DLC0243"])
counttry <- (counttry[,colnames(counttry)!="DLC0219"])


counttry<- countsTable[vec<3]
vec=vec[vec<3]

nsamp= length(counttry)
vec <- vector(mode="integer", length=length(counttry[1:nsamp]))
for (i in 1:length(counttry[1:nsamp])){
  vec[i]=categories$V11[which(categories$V1==colnames(counttry)[i])]
}


