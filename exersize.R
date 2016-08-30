library('RTCGAToolbox')
load('ann_chr3.RData')

getFirehoseDatasets()
getFirehoseRunningDates(last = NULL)
readData = getFirehoseData (dataset="UVM", runDate="20151101",forceDownload = TRUE,
                            Clinic=FALSE, Mutation=FALSE, Methylation=FALSE, 
                            RNAseq2_Gene_Norm=TRUE)

mRNA <- getData(readData, 'RNASeq2GeneNorm')

gs <- read.delim('custom.bindea.gmt')

######################################################
##http://stackoverflow.com/questions/6602881/text-file-to-list-in-r

# Read in the data
x <- scan("custom.bindea.gmt", what="", sep="\n")
# Separate elements by tab
y <- strsplit(x, "\t")
# Extract the first vector element and set it as the list element name
names(y) <- sapply(y, `[[`, 2)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
y <- lapply(y, `[`, -1:-2)

#y <- lapply(y, function(x) x[-1]) # same as above
es <- gsva(mRNA, y, method = 'ssgsea', rnaseq = TRUE )
es_transpose_scale <- scale(t(es), center = TRUE, scale = TRUE)
rownames(es_transpose_scale) <- gsub('-', '\\.', gsub('-...-...-....-..', '' , rownames(es_transpose_scale)))

ann <- ann[match(rownames(es_transpose_scale), rownames(ann)),]


library(NMF)
aheatmap(es_transpose_scale,
         hclustfun=function(d) hclust(d, method="ward.D2"),
         annRow = ann)

library(Cairo)
CairoPDF(file = 'cibersort_TCGA.pdf',
         width =7.5, height = 7.5, pointsize = 16)
#layout(matrix(c(1,1,2,2), ncol = 2, byrow = TRUE),
#       widths = c(1,1),
#       heights = c(276,95)) 
#########################################################################

ann <- data.frame(Chr3 = Chr3)

ah<- aheatmap(data_ciber[,252:273], 
              hclustfun=function(d) hclust(d, method="complete"),
              #Colv = colv,
              annRow = ann#,
              #annCol = ann_col,
              #annColors = ann_colors_e,
              #cex = 2,
              #labRow = rep('',dim(data_hc_e)[1]),
              #labCol = rep('',dim(data_hc_e)[2])
              #fontsize = 12,
              #cexCol = ,
              #fontsize = 16,
              #labCol = rep('',dim(data_hc_e)[2]),
              #Colv = colv
              #reorderfun = function(d, w) reorder(d, 10)
)
dev.off()
mRNA[1:10, 1:10]

sum(duplicated(gsub('-...-...-....-..', '', colnames(mRNA))))
colnames(mRNA) <- gsub('-...-...-....-..', '', colnames(mRNA))

write.table(mRNA, file = 'TCGA_mixture.txt', append = FALSE, 
            quote = FALSE, sep = '\t',
            row.names = TRUE,
            col.names = TRUE 
)
