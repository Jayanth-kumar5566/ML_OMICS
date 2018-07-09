library(tgp)
library(caret)
library(corrplot)

data=read.csv("Documents/MS/ML_OMICS/Clustering/Luminex/Luminex.csv")
data=na.omit(data)
rownames(data)=data$PatientID
data=data[,2:43]

trainIndex <- createDataPartition(data$BSI_Bronchiectasis_severity_index, p = 0.8,list = FALSE,times = 1)


X_train <- data[trainIndex,c(1:41)]# X
Z_train <- data[trainIndex,42]

X_test <-data[-trainIndex,c(1:41)]
Z_test <-data[-trainIndex,42]

# # use a Treed Linear Model (see documentation for more information)
#fit <- btgpllm(X=X_train, Z=Z_train, XX=X_test)
fit <- btgpllm(X=X_train, Z=Z_train, XX=X_test)#,basemax=47)
errors <- sqrt(mean((fit$ZZ.mean - Z_test)^2))
tgp.trees(fit)
errors


corrplot(cor(data))


#Clustering
library(cluster)

#kMeans
wit=c(NA)
for(i in 2:10){wit[i]=kmeans(data,i)$tot.withinss}
plot(1:10,wit)
lines(1:10,wit)

#hclust
hc=hclust(dist(data),method="ward.D2")
plot(hc)
