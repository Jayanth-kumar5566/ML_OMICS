#Importing Modules
library("SNFtool")
library("vegan")
library("igraph")
library("Rtsne")
#Importing the data
b_data=read.csv("Documents/MS/ML_OMICS/SNF_standard/bacteria.csv")
f_data=read.csv("Documents/MS/ML_OMICS/SNF_standard/fungi.csv")
#Bray Curtis Similarity
b_dsim=vegdist(b_data[,-1],method='bray',diag=TRUE,upper=TRUE)
f_dsim=vegdist(f_data[,-1],method='bray',diag=TRUE,upper=TRUE)
W1=(as.matrix(b_dsim)-1)*-1
W2=(as.matrix(f_dsim)-1)*-1
#plots graphs
x=W2
ig <- graph.adjacency(x, mode="undirected", weighted=TRUE)
labels=spectralClustering(x,2)
V(ig)$lab=labels
colrs <- c("red", "blue")
V(ig)$color <- colrs[V(ig)$lab]
E(ig)$width <- E(ig)$weight*0.1
ig <- delete_edges(ig, E(ig)[weight>0.99999])
rbPal <- colorRampPalette(c('white','grey','green','orange','blue'))
E(ig)$color<-rbPal(10)[as.numeric(cut(E(ig)$weight,breaks = 10))]
#ig <- delete_edges(ig, E(ig)[weight<0.85])
la=Rtsne(as.matrix(b_dsim),pca=FALSE,is_distance = TRUE)
lay=la$Y

lay=cbind(labels*10,labels*10)
for (i in 1:229){for (j in 1:2){lay[i,j]=lay[i,j]+runif(1,0,7)}}

plot.igraph(ig,vertex.label=NA,vertex.size=5,layout=lay)
legend("topleft",title="Fungal",legend=c(1:10),col =rbPal(10),pch=20)
#SNF
W = SNF(list(W1,W2),23,15)
#Spectral Clustering
estimateNumberOfClustersGivenGraph(W)
displayClusters(W, spectralClustering(W, K = 2))
#Results
C=2
labels=spectralClustering(W,C)
#Visuvalization
h=max(W[lower.tri(W)])
l=min(W[lower.tri(W)])
W=(W-l)/(h-l)
diag(W)<-1.2
ig <- graph.adjacency(W, mode="undirected", weighted=TRUE)
V(ig)$lab=labels
colrs <- c("red", "blue")
V(ig)$color <- colrs[V(ig)$lab]
E(ig)$width <- E(ig)$weight*2
ig <- delete_edges(ig, E(ig)[weight>1.1])
#ig <- delete_edges(ig, E(ig)[weight<0.85])
rbPal <- colorRampPalette(c('white','grey','green','orange','blue'))
E(ig)$color<-rbPal(10)[as.numeric(cut(E(ig)$weight,breaks = 10))]
la=Rtsne(1.2-W,pca=FALSE,is_distance = TRUE)
lay=la$Y

lay=cbind(labels*10,labels*10)
for (i in 1:229){for (j in 1:2){lay[i,j]=lay[i,j]+runif(1,0,7)}}

plot.igraph(ig,vertex.label=NA,vertex.size=5,layout=lay)
legend("topleft",title="SNF",legend=c(1:10),col =rbPal(10),pch=20)
#Writing to a file
write.csv(W,"Documents/MS/ML_OMICS/SNF_standard/merged_matrix.csv")
write.csv(labels,"Documents/MS/ML_OMICS/SNF_standard/labels.csv")

