########################################
## Pancancer Immune ####################

library('ssGSEA4TCGA')

rdate <- getFirehoseRunningDates(last = NULL)

dset <- getFirehoseDatasets()
dset <- dset[c(2,3,6,12,13,21,22,24,30)]#,31,34)]
dset <- c('BRCA', 'LUAD')
gs <- gs_gmt('custom.bindea.gmt')

data <- pancancer_ssGSEA(dset)
data_t <- scale(t(data[[1]]), center = TRUE, scale = TRUE)

TvsN <- c(rep('T', dim(data_t)[1]))
TvsN[grep('^TCGA-..-....-1.*', rownames(data_t))] <- 'N'
TvsN <- factor(TvsN, levels = c('T', 'N'),
               labels = c('Tumor', 'Normal'))
summary(TvsN)

ann <- data.frame(class = as.factor(results_row[[6]]$consensusClass),
                  type = data[[2]], 
                  Tumor = TvsN)

library(ConsensusClusterPlus)

results_col = ConsensusClusterPlus(data_t,maxK=5,reps=5000,pItem=0.8,pFeature=1,
                                   title='consensus_col',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

results_row = ConsensusClusterPlus(data[[1]],maxK=10,reps=500,pItem=0.8,pFeature=1,
                                   title='consensus_by_row',
                                   clusterAlg="hc",
                                   innerLinkage = "ward.D2",
                                   finalLinkage = "ward.D2",
                                   distance="euclidean",
                                   plot="pdf")

library(NMF)
library(Cairo)
ann_col <- data.frame(immune_class = factor(results_col[[3]]$consensusClass))

CairoPDF(file = 'pancancer_TCGA.pdf',
         width =7.5, height = 7.5, pointsize = 16)
aheatmap(data_t,
         hclustfun=function(d) hclust(d, method="ward.D2"),
         annRow = ann,
         annCol = ann_col,
         Colv = results_col[[3]]$consensusTree,
         Rowv = results_row[[6]]$consensusTree,
         labRow = rep('',dim(data_t)[1]))
dev.off()
