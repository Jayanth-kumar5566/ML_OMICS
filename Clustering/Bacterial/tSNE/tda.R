library(vegan)
library(Rtsne)
library(cluster)

data=read.csv("Documents/MS/ML_OMICS/Clustering/Bacterial/bacterial_microbiome.csv")
result <- data[-1]
row.names(result) <- data$id
#Remove rows with missing values
result<-na.omit(result)
#Work with result
Bray_curtis_diss=vegdist(result,method="bray",na.rm = TRUE)

Diss_matrx<-as.matrix(Bray_curtis_diss)
tsne_obj <- Rtsne(Diss_matrx,dims=2,perplexity=30,pca=FALSE,theta=0.0,max_iter=2000,
                  verbose=TRUE,is_distance = TRUE)
plot(tsne_obj$Y)

#----------------PAM Analysis----------------
sil_width <- c(NA)

for(i in 2:15){
  pam_fit <- pam(Bray_curtis_diss,
                 diss = TRUE,
                 k = i)
  sil_width[i] <- pam_fit$silinfo$avg.width
}

plot(1:15, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:15, sil_width)

#the silhouette anlaysis yeilds that there are 3 clusters(from plot)

pam_fit <- pam(Bray_curtis_diss, diss = TRUE, k = 5)

#Visuvalization on tSNE
plot(tsne_obj$Y,col=pam_fit$clustering)

#Silhoutee width of clusters
pam_fit$silinfo$clus.avg.widths

#Some statistics of the clusters
pam_fit$clusinfo

#----------Hirechial Clustering------------------------
hc.m = hclust(Bray_curtis_diss, method="median")
hc.s = hclust(Bray_curtis_diss, method="single")
hc.c = hclust(Bray_curtis_diss, method="complete")
windows(height=3.5, width=9)
layout(matrix(1:3, nrow=1))
plot(hc.m)
plot(hc.s)
plot(hc.c)
#----------------DBSCAN---------------------------
layout(matrix(1:2, nrow=1))
plot(density(na.omit(Bray_curtis_diss[upper.tri(Bray_curtis_diss)])), main="kernel density")
plot(ecdf(Bray_curtis_diss[upper.tri(Bray_curtis_diss)]), main="ECDF")

