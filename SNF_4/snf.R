library("SNFtool")
library("vegan")
b_data=read.csv("Documents/MS/ML_OMICS/SNF/bacteria.csv")
f_data=read.csv("Documents/MS/ML_OMICS/SNF/fungi.csv")
b_dsim=vegdist(b_data[,-1],method='bray',diag=TRUE,upper=TRUE)
f_dsim=vegdist(f_data[,-1],method='bray',diag=TRUE,upper=TRUE)
W1=(as.matrix(b_dsim)-1)*-1
W2=(as.matrix(f_dsim)-1)*-1
displayClusters(W1, spectralClustering(W1, K = 4))
W = SNF(list(W1,W2),20,100)
estimateNumberOfClustersGivenGraph(W)
displayClusters(W, spectralClustering(W, K = 2))
C=4
labels=spectralClustering(W,C)
write.csv(W,"Documents/MS/ML_OMICS/SNF_4/merged_matrix.csv")
write.csv(labels,"Documents/MS/ML_OMICS/SNF_4/labels.csv")

for (i in 1:229){
  W = SNF(list(W1,W2),i,100)
  z=estimateNumberOfClustersGivenGraph(W)$`Eigen-gap best`
  print(paste(c(i,z),sep = " "))
  labels=spectralClustering(W,z)
  write.csv(labels,paste(c("/tmp/LAB/labels",i,".csv"),collapse = ""))
  write.csv(W,paste(c("/tmp/LAB/matrix",i,".csv"),collapse = ""))
}
