########################################
## Pancancer Immune ####################

library('ssGSEA4TCGA')

rdate <- getFirehoseRunningDates(last = NULL)

dset <- getFirehoseDatasets()
dset <- dset[c(2,3,6,12,13,21,22,24,30)]#,31,34)]

gs <- gs_gmt('custom.bindea.gmt')

data <- pancancer_ssGSEA(dset)
data_t <- scale(t(data[[1]]), center = TRUE, scale = TRUE)

TvsN <- c(rep('T', dim(data_t)[1]))
TvsN[grep('^TCGA-..-....-1.*', rownames(data_t))] <- 'N'
TvsN <- factor(TvsN, levels = c('T', 'N'),
               labels = c('Tumor', 'Normal'))
summary(TvsN)

ann <- data.frame(type = data[[2]], Tumor = TvsN)

library(NMF)
aheatmap(data_t,
         hclustfun=function(d) hclust(d, method="ward.D2"),
         annRow = ann)

