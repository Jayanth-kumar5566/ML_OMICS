#Functions to source into the app
library(SNFtool)
library("vegan")
library("reticulate")
use_python("/usr/local/bin/python")
source_python("Python_codes/plot_mat.py")
source_python("Python_codes/sil.py")
#=================Trail data====================
data<-read.csv("./../Data/bacteria.csv",header = TRUE,row.names = 1)
#=====Function for plotting individual biomes===================
biome_plot<-function(data,k){
  dsim=vegdist(data,method='bray',diag=TRUE,upper=TRUE)
  sim=(as.matrix(dsim)-1)*-1
  labels=spectralClustering(sim,k)
  labels=as.data.frame(labels,row.names = row.names(data))
  sim=as.data.frame(sim)
  #bplot(sim,labels)
  m<-bplot(sim,labels)
  return(heatmap(as.matrix(m),Rowv = NA, Colv = NA, scale = "none",
                 main = "Spectral clustering",xlab = "Sample ID",ylab = "Sample ID"))
}

#==============Function for giving silhouette score================
#Given Similarity matrix and labels
biome_silhouette<-function(sim,labels){
  silhouette_score(sim,labels)
}
#==================================================================