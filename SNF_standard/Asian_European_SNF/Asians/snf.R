library("SNFtool")
library("vegan")
#library("RColorBrewer")
b_data=read.csv("bacteria.csv")
f_data=read.csv("fungi.csv")
#l_data=read.csv("Documents/MS/ML_OMICS/SNF_standard/luminex.csv")
b_dsim=vegdist(b_data[,-1],method='bray',diag=TRUE,upper=TRUE)
f_dsim=vegdist(f_data[,-1],method='bray',diag=TRUE,upper=TRUE)
#l_dsim=vegdist(l_data[,-1],method='bray',diag=TRUE,upper=TRUE)
W1=(as.matrix(b_dsim)-1)*-1
W2=(as.matrix(f_dsim)-1)*-1
#W3=(as.matrix(l_dsim)-1)*-1
x=W2
estimateNumberOfClustersGivenGraph(x)
displayClusters(x, spectralClustering(x, K = 2))
W = SNF(list(W1,W2),13,15)

#Divide W by 0.5 W/0.5

estimateNumberOfClustersGivenGraph(W)
displayClusters(W, spectralClustering(W, K = 2))
C=2
#W=x
labels=spectralClustering(W,C)
write.csv(W,"merged_matrix.csv")
write.csv(labels,"labels.csv")

for (i in 1:229){
  W = SNF(list(W1,W2),i,100)
  z=estimateNumberOfClustersGivenGraph(W)$`Eigen-gap best`
  print(paste(c(i,z),sep = " "))
  labels=spectralClustering(W,z)
  write.csv(labels,paste(c("/tmp/LAB/labels",i,".csv"),collapse = ""))
  write.csv(W,paste(c("/tmp/LAB/matrix",i,".csv"),collapse = ""))
}


#-----------For Luminex----
l_data=read.csv("luminex.csv")
l_diss=dist(l_data[,-1],method="euclidean",diag=TRUE,upper=TRUE)
#Dist1 = (dist2(as.matrix(l_data[,-1]),as.matrix(l_data[,-1])))^(1/2)
#l_sim1=affinityMatrix(Dist1,K=22,sigma=0.5)
l_sim=affinityMatrix(as.matrix(l_diss),K=22,sigma=0.5)
W=l_sim
estimateNumberOfClustersGivenGraph(W)
displayClusters(W, spectralClustering(W, K = 2))
C=2
labels=spectralClustering(W,C)
write.csv(W,"merged_matrix.csv")
write.csv(labels,"labels.csv")
