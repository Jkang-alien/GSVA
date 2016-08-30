########################################
## Pancancer Immune ####################

library('RTCGAToolbox')
library(GSVA)

dset <- getFirehoseDatasets()
dset <- dset[c(2,3,6,12,13,21,22,24,30)]#,31,34)]
rdate <- getFirehoseRunningDates(last = NULL)

######################################################
##http://stackoverflow.com/questions/6602881/text-file-to-list-in-r

gs_gmt <- function(a){
# Read in the data
x <- scan(a, what="", sep="\n")
# Separate elements by tab
y <- strsplit(x, "\t")
# Extract the first vector element and set it as the list element name
names(y) <- sapply(y, `[[`, 2)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
y <- lapply(y, `[`, -1:-2)
return(y)
}

ssGSEA <- function (a){
  readData = getFirehoseData (dataset=a, runDate=rdate[1], RNAseq2_Gene_Norm=TRUE)
  mRNA <- getData(readData, 'RNASeq2GeneNorm')
  es <- gsva(mRNA, gs, method = 'ssgsea', rnaseq = TRUE )
  return (es)
}


pancancer_ssGSEA <- function(a){
  d <- c()
  for (i in a){
    b <- ssGSEA(i)
    d <- cbind(d, b)
  }
  return (d)
}

gs <- gs_gmt('custom.bindea.gmt')
data <- pancancer_ssGSEA(dset)
site <- c(rep('colorectal',263), rep('luad', 576))
ann <- data.frame(site = site)

  
data_t <- scale(t(data), center = TRUE, scale = TRUE)
library(NMF)
aheatmap(data_t,
         hclustfun=function(d) hclust(d, method="ward.D2"),
         annRow = ann)
