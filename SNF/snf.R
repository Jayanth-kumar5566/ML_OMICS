library("SNFtool")
library("vegan")
b_data=read.csv("Documents/MS/ML_OMICS/SNF/bacteria.csv")
f_data=read.csv("Documents/MS/ML_OMICS/SNF/fungi.csv")
b_dsim=vegdist(b_data[,-1],method='bray',diag=TRUE,upper=TRUE)
f_dsim=vegdist(f_data[,-1],method='bray',diag=TRUE,upper=TRUE)
W1=(as.matrix(b_dsim)-1)*-1
W2=(as.matrix(f_dsim)-1)*-1
estimateNumberOfClustersGivenGraph(W)
displayClustersWithHeatmap(W1, spectralClustering(W1, K = 4))
W = SNF(list(W1,W2),20,100)
displayClustersWithHeatmap(W, spectralClustering(W, K = 2))
C=2
labels=spectralClustering(W,C)
