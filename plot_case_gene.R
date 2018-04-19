##### Author: S Saberi
##### Version: 1.1.0
##### Last modified: Aug 30, 2017
##### Function: multi hits impact case gene plot ###########
##### Input: 3-column data frame with header: case, gene, type or a directory that contains these files
##### Output: case-gene-aberration plot (oncoplot)
##### Usage: Rscript plot_case_gene.R test_inFile.txt NULL test_plot.pdf
##### Test 4 column file: Rscript plot_case_gene.R ./test_geneSampleOrder2.txt test_newplot.pdf ./test_inFile2.txt testout.txt 
##### Test 5 column file (includes frequency): Rscript plot_case_gene.R NULL test_geneSampleOrder3.txt test_oncoplot.pdf ./test_inFile3.txt testout.txt
##### Test 5 column file (includes frequency): Rscript plot_case_gene.R test_geneSampleOrder3.txt test_geneSampleOrder3.txt test_oncoplot.pdf ./test_inFile3.txt testout.txt
##### Test 5 column file (includes frequency): Rscript plot_case_gene.R NULL NULL test_oncoplot.pdf ./test_inFile3.txt testout.txt

add.alpha <- function(col, alpha = 1){
  # Add an alpha value to a colour
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}



PlotOncoplotPerPage <- function (case.gene.type, plot.name, genes.per.page = 1000, sample.order, sample.color, genes = NULL) {
  # Defines the order of the genes in the oncoplot and plots them page by page.
  #
  # Args:
  #   case.gene.type: 3-column data frame with header: case, gene, type.
  #   genes.per.page: number of genes (rows in oncoplot) plotted per page. 
  #   sample.order: a list: $sample_order - a unique list of samples for labeling each column of the plot 
  #                         $color - a unique list of color for samples according to their type - same length as sample_order
  #
  # Returns:
  #   NULL.
  
  # Error handling
  if (nrow(case.gene.type) <=1) {
    stop(" case.gene.type has no rows.")
  } else if (ncol(case.gene.type) == 3) {
    if (colnames(case.gene.type) != c("case", "gene", "type")) {
      stop(" case.gene.type has 3 columns but missing either of the case, gene or type columns.")
    } 
  } else if(ncol(case.gene.type) == 4) {
    if (colnames(case.gene.type) != c("case", "gene", "type", "project")) {
      stop(" case.gene.type has 4 columns but missing either of the case, gene, type or project columns.")
    } 
  } else if(ncol(case.gene.type) == 5) {
    if (colnames(case.gene.type) != c("case", "gene", "type", "project", "frequency")) {
      stop(" case.gene.type has 5 columns but missing either of the case, gene, type, project or frequency columns.")
    } 
  }
  
  case.gene <- names(sort(table(unique(paste(case.gene.type[,"case"], case.gene.type[, "gene"], sep = ":"))), 
                          decreasing = TRUE))
  
  getGene <- function(x) {gname <- strsplit(x, split = ":")[[1]][2];return(gname)}
  tbl <- sort(sort(table(sapply(case.gene, getGene))), decreasing = TRUE)
  genes.uniq <- names(tbl)
  Freq <- as.numeric(tbl)


  tbl <- data.frame(names = genes.uniq, Freq = Freq)
  tbl <- tbl[order(Freq),]
  print('aha')
  print(tbl)
  pdf('C_T-Freq.oncoplot.pdf')
    p <- ggplot(data=tbl, aes(x=reorder(names,-Freq), y=Freq) ) + geom_bar(stat = "identity",fill= c("blue")) +
    theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
    p <- p + ylab('Mutation freq %') + xlab('')
    print(p)
  dev.off()
  n.genes <- length(genes.uniq)
  pdf(plot.name, width=20, height=10)
  # Iterate through batches of genes and plot the oncoplot
  for (k in 1:(floor(n.genes/genes.per.page)+1)) {  
    # Define gene order
    
    gene.order <- genes.uniq[((k - 1) * genes.per.page + 1) : (k * genes.per.page)]
    # Define sample order
    
    if (is.null(sample.order)) { # data-driven sample order per page
      sample_order <- c()
      
      for (g in gene.order) {
        sample_order <- c(sample_order, case.gene.type[case.gene.type$gene == g, "case"])
      }
      sample_order <- union(sample_order, case.gene.type[, "case"])
      sample_order <- sample_order[!duplicated(sample_order)]
      
    } else { # user-specified sample order/color

      
      color = as.character(sample.order$color)
      sample_order.tmp=c()
      for (c in unique(sample.order$color)) { 
        samples <- sample.order[sample.order$color==c,]$sample_order
        for (g in gene.order[!is.na(gene.order)]) {
          case.gene.type1 <- case.gene.type[case.gene.type$gene == g, ]
          case.gene.type1 <- intersect(case.gene.type1$case , samples)
          
          sample_order.tmp <- c(sample_order.tmp, case.gene.type1 )
        }
        sample_order.tmp <- union(sample_order.tmp, sample.order[sample.order$color==c,]$sample_order )
        sample_order.tmp <- sample_order.tmp[!duplicated(sample_order.tmp)]
      }
      #sample_order.tmp <- union(sample_order.tmp, case.gene.type[, "case"])
      sample_order.tmp <- sample_order.tmp[!duplicated(sample_order.tmp)] 
      sample_order  <-  as.character(sample_order.tmp)
    }
    
    # Define sample color
    if (is.null(sample.color)) {
      color = rep("blue", length(sample_order))
    } else {
      rownames(sample.color) <- sample.color$sample_order
      color = as.character(sample.color[sample_order, "color"])
      
    }
    if(!is.null(genes)) {gene.order <-genes }
    gene.sample.order = list(gene_order = gene.order, sample_order = sample_order, color = color)
    
    # Plot the oncoplot for each page
    ggplot.obj <- PlotOncoplot(case.gene.type, gene.sample.order, genes = genes)
    print(ggplot.obj)
  }
  dev.off()
  return(NULL)
}







PlotOncoplot <- function (case.gene.type, gene.sample.order = NULL, genes = NULL) {
  # Plots an oncoplot based on the specific sample order (if any).
  #
  # Args:
  #   case.gene.type: 3-column data frame with header: case, gene, type.
  #   sample.order: a list:$gene_order - order of genes in the oncoplot  
  #                        $sample_order - a unique list of samples for labeling each column of the plot
  #	                       $color - a unique list of color for samples according to their type - same length as sample_order
  #
  # Returns:
  #   ggplot.obj: ggplot object
  
  # Define gene_order, sample_order and sample label color if it is not assigned
  if (is.null(gene.sample.order)) {     
    gene_order = unique(case.gene.type[ ,"gene"])
    sample_order = unique(case.gene.type[ ,"case"])
    color = rep("blue", length(sample_order))
    gene.sample.order = list(gene_order = gene_order, sample_order = sample_order, color = color)
  }
  
  # Assign color to each type --------------------------
  get_palette <- colorRampPalette(brewer.pal(9, "Set1"))
  types <- unique(case.gene.type[ ,"type"]) ## types = c("type1","type2")
  type.cols <- get_palette(length(types))   # Set the main colors for each aberration type
  names(type.cols) <- types  
  
  if (ncol(case.gene.type) == 5) {  # Set alpha for colors if the frequency column is present 
    case.gene.type[ ,"type"] <- paste(case.gene.type[ ,"type"], case.gene.type[ ,"frequency"], sep = ":") # New type column
    types.freq <- unique(case.gene.type[ ,"type"]) # Define new type (combination of type and frequency) ## types.freq = c("type1_2","type1_5","type2_1","type2_3","type2_4","type2_5")
    
    type.cols.alpha = rep(0, length(types.freq))  # New color scheme (alpha scaled)
    names(type.cols.alpha) <- types.freq

    
    
    for (type in types) {  # Iterate through main aberration types and alpha scale based on frequency 
      type.freq.tmp <- grep(type, types.freq, value = TRUE)
      for (ty.t in type.freq.tmp) {
        color.tmp <- type.cols[strsplit(ty.t, split = ":")[[1]][1]]
        alpha.tmp <- as.numeric(strsplit(ty.t, split = ":")[[1]][2])
        if (alpha.tmp > 5) { alpha.tmp = 5} # Maximum alpha value
        type.cols.alpha[ty.t] <- add.alpha(color.tmp, alpha = alpha.tmp/5)
      }
    }
    
    type.cols <- type.cols.alpha  # Use alpha scaled colors
  }

  
  case.gene.type <- case.gene.type[case.gene.type[,"gene"]%in%gene.sample.order$gene_order,]

  # Reformat boxes for multi hits ---------------------------
  dt <- data.table(unique(case.gene.type))
  dt[, var1 := (match(gene, gene.sample.order$gene_order))]
  dt[, var2 := match(case, gene.sample.order$sample_order)]
  dt[, shift := (1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(var1, var2)]
  dt[, height := 1/.N, by = list(var1, var2)]
  dt$type = factor(dt$type, levels = unique(dt$type))


  #write.table(dt, file=data_outfile, append=TRUE, row.names=FALSE, quote=FALSE)  
  
  ## ### define sample label color	---------------------------
  ## sample_cols<-unique(data.frame(sample=dt$case, color=dt$color))
  ## sample.cols<-as.vector(sample_cols[,"color"])
  ## names(sample.cols)<-as.vector(sample_cols[,"sample"])

  for (case in unique(as.character(gene.sample.order$sample_order))) {
    if (length(intersect(case, dt$case))==0){
      print (case)
      tmp <- dt[length(dt$case),]
      tmp$case <- case
      tmp$var2 <- NA
      dt <- rbind(dt,tmp)
      print(tmp)

    }
  }
  
  ggplot.obj <- ggplot(dt, aes(var2, y = var1 + shift, fill = type, height = height)) + geom_tile(color = "gray95", size = 1.) + 
    scale_fill_manual(values = type.cols) + xlab("Case ID") + ylab("Gene") + 
    scale_x_continuous(breaks = seq(1:length(gene.sample.order$sample_order)+1), labels = gene.sample.order$sample_order) + 
    scale_y_reverse(breaks=seq(1:length(gene.sample.order$gene_order)), labels = gene.sample.order$gene_order)  + 
    #coord_cartesian(xlim = seq(1:length(gene.sample.order$sample_order)), labels = gene.sample.order$sample_order) + 
    theme_bw() + 
    theme(text = element_text(size=16), axis.text.x = element_text(angle = 90, hjust = .5, color = gene.sample.order$color) )  + 
    coord_fixed(ratio=1) 
  
  return(ggplot.obj)
}

##### 
##### MAIN BODY ---------------------------
##### 
library("RColorBrewer")
library("data.table")
library("ggplot2")

# Handling arguments
# options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
#print(args)

sample.order.file<- args[1]
sample.color.file<- args[2]
plot.name <- args[3]
snv_infile <- args[4] 
genes <- args[5] 

print(args)
#data_outfile <- args[5]


#if (file.exists(data_outfile)) file.remove(data_outfile)

if (sample.order.file == "NULL") { # read sample order from a file
  sample.order = NULL
} else {
  sample.order <- read.delim(sample.order.file)
}
if (sample.color.file == "NULL") { # read sample color from a file
  sample.color = NULL
} else {
  sample.color <- read.delim(sample.color.file)
}
if (genes == "NULL") { # read sample color from a file
  genes = NULL
} else {
  genes <- read.table(genes,header = T,as.is = T,sep = '\t')
  genes <- genes$gene
}

#load_snv_data
case.gene.type = c()
data.snv <- read.table(snv_infile, as.is=T, header=T, check.names = FALSE, sep="\t")
case.gene.type <- rbind(case.gene.type, data.snv)
remove(data.snv)


PlotOncoplotPerPage(case.gene.type, paste(plot.name, '.pdf', sep=""), genes.per.page = 200, sample.order = sample.order, 
                    sample.color = sample.color, genes)




