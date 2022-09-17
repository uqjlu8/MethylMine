#Mixomics data integration and analysis

#load library
library(mixOmics)

#load data
data(breast.TCGA)
# extract training data and name each data frame
X <- list(mRNA = breast.TCGA$data.train$mrna, 
          miRNA = breast.TCGA$data.train$mirna, 
          protein = breast.TCGA$data.train$protein)
Y <- breast.TCGA$data.train$subtype

#check summary
summary(Y)

#generate list
list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), protein = c(5, 5))

#plot omic-specific data
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
plotIndiv(MyResult.diablo) ## sample plot

#correlation circos plots
plotVar(MyResult.diablo) ## variable plot

plotIndiv(MyResult.diablo, 
          ind.names = FALSE, 
          legend=TRUE, cex=c(1,2,3),
          title = 'DIABLO - TCGA-BRCA')
          
          
#Correlation plot 
plotDiablo(MyResult.diablo, ncomp = 1)


#bloack plsda
MyResult.diablo2 <- block.plsda(X, Y)

#heatmap
# minimal example with margins improved:
# cimDiablo(MyResult.diablo, margin=c(8,20))
# extended example:
cimDiablo(MyResult.diablo, color.blocks = c('darkorchid', 'brown1', 'lightgreen'), comp = 1, margin=c(8,20), legend.position = "right")

#plotLoadings(MyResult.diablo, contrib = "max")
plotLoadings(MyResult.diablo, comp = 2, contrib = "max")
