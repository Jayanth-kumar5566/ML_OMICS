library(vegan)
library(SNFtool)
library(caret)

#label=c() #import [must be in numbers like 1,2,3]

file_n=c("Data/SMOTE/data_25.csv","Data/SMOTE/data_50.csv","Data/SMOTE/data_75.csv","Data/SMOTE/data_100.csv","Data/SMOTE/data_smoted.csv")

F_Acc=c()

for (j in file_n){
print(j)
data=read.csv(j)

b_data=data[,1:134]
f_data=data[,135:539]

label=data$exacgrp
label=ifelse(label=="yes",2,1)
Acc=c()

for (i in 1:100){
  set.seed(runif(1,0,1e9))
  data_set_size <- floor(nrow(b_data)/3)
  indexes <- sample(1:nrow(b_data), size = data_set_size)
  b_data_r<-rbind(b_data[-indexes,],b_data[indexes,])
  f_data_r<-rbind(f_data[-indexes,],f_data[indexes,])
  train_groups = label[-indexes]
  test_groups = label[indexes]
  
  #Bray Curtis Similarity
  b_dsim=vegdist(b_data_r[,-1],method='bray',diag=TRUE,upper=TRUE)
  f_dsim=vegdist(f_data_r[,-1],method='bray',diag=TRUE,upper=TRUE)
  W1=(as.matrix(b_dsim)-1)*-1
  W2=(as.matrix(f_dsim)-1)*-1
  
  dataL=c()
  dataL[[1]]=W1
  dataL[[2]]=W2
  
  
  grpPredict<- function (dataL, groups, K = 22, t = 15, method = 1, type="pred") 
  {
    #Assumes that dataL=rbind(train,test), ie similarity matrix is calculated in this way
    #Assumes the dataL are bray-curtis similarities
    W <- SNF(dataL , K, t)
    Y0 <- matrix(0, nrow(dataL[[2]]), max(groups)) #2 or 1 doesn't matter
    for (i in 1:length(groups)) {
      Y0[i, groups[i]] <- 1
    }
    Y <- SNFtool:::.csPrediction(W, Y0, method)
    newgroups <- rep(0, nrow(dataL[[1]])) #2 or 1 doesn't matter
    if (type=="prob"){
      return(Y)
    }
    else{
      for (i in 1:nrow(Y)) {
        newgroups[i] <- which(Y[i, ] == max(Y[i, ]))
      }
      return(newgroups)
    }
  }
  
  #res<-grpPredict(dataL,train_groups,type="prob")
  prediction_for_table<-grpPredict(dataL,train_groups,type="pred")
  x<-tail(prediction_for_table,length(test_groups))
  print(x)
  cs<-confusionMatrix(as.factor(x),as.factor(test_groups))
  Acc=c(Acc,cs$overall[1])
}

print(mean(Acc))

F_Acc=c(F_Acc,mean(Acc))
}

#Class Balancing 
#indexes of labels which are 1 then sample 65 from them and combine them

plot(c(25,50,75,100,152),F_Acc)
