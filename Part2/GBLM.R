#Reading the dataset
data=read.csv("data_combined.csv")
rownames(data)<-data$PatientID
data$PatientID<-NULL

c1_data=subset(data,x==1)
c2_data=subset(data,x==2)
c1_data$x<-NULL
c2_data$x<-NULL

#Checking for relative abundance structure
r_s=c()
for (i in 1:229){
  r_s=c(r_s,sum(data[i,]))
}
print(var(r_s))
remove(r_s,i,data)

library(mboost)
cr_mat=cor(c1_data,method = "spearman")
cr_mat[is.na(cr_mat)]<-0

cols=which(abs(cr_mat[1,]) > 0.05, arr.ind = T)
rownames(as.data.frame(cols))
