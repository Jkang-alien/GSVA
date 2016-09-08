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

clinic_luad <- getFirehoseData('LUAD', rdate[1], Clinic = TRUE)
clinic_luad_data <- getData(clinic_luad, 'Clinical')

clinic_luad_data <- survivalTCGA(clinic_luad_data)

library(survival)
library(rms)


diff = survdiff(Surv(surv_months, vital_status == 1)~ case_class, 
                data = data)
diff

svg(file = "Figure3.svg", pointsize = 10,
    width = 7.5 , height = 4,)
layout(matrix(c(1,2), ncol = 2, byrow = TRUE))
par(mar=c(5,3,1,4), mgp = c(2, 1, 0))

fit = npsurv(Surv(surv_months, vital_status == 1)~ case_class, 
             data = data)

strata = levels(data$case_class)

survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:6),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=data$case_class,
         col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.3, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:6), col=data$case_class, cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*10, 0.8, 'P-value: 0.052', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')


cox <- coxph(Surv(surv_months, vital_status == 1)~ case_class, 
             data = data)

summary(cox)


