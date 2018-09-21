library(vegan)
library(SNFtool)
library(caret)

#label=c() #import [must be in numbers like 1,2,3]

b_data=read.csv("SMOTE/bacteria_smoted.csv")
f_data=read.csv("SMOTE/fungi_smoted.csv")
l=read.csv("SMOTE/DM_CM_smoted.csv")
label=c(l$exacgrp)

Acc=c()

for (i in 1:100){
set.seed(runif(1,0,1e9))

#Under-sampling
# table(label)
# ind=which(label==1)
# ind_sam=sample(ind,65)
# ind_low=which(label==2)
# ind_tot=c(ind_sam,ind_low)
# 
# b_data=b_data[ind_tot,]
# f_data=f_data[ind_tot,]
# label=label[ind_tot]
  

data_set_size <- floor(nrow(b_data)/3)
indexes <- sample(1:nrow(b_data), size = data_set_size)
b_data_r<-rbind(b_data[-indexes,],b_data[indexes,])
f_data_r<-rbind(f_data[-indexes,],f_data[indexes,])
train_groups = label[-indexes]
test_groups = label[indexes]
print(table(train_groups))

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
  #W <- dataL[[2]] #1 or 2 depending the bact or fung
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

res<-grpPredict(dataL,train_groups,type="prob")
prediction_for_table<-grpPredict(dataL,train_groups,type="pred")
x<-tail(prediction_for_table,length(test_groups))
print(x)
cs<-confusionMatrix(as.factor(x),as.factor(test_groups))
Acc=c(Acc,cs$overall[1])
}

print(mean(Acc))
#Class Balancing 
#indexes of labels which are 1 then sample 65 from them and combine them