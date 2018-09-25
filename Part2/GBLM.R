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
# For c1_data
cr_mat=cor(c1_data,method = "spearman") #Compute Spearman Correlation
cr_mat[is.na(cr_mat)]<-0  # Assign 0 to NA values(NA due to zero STD)

i=1
#for (i in 1:539){
  x_nam=rownames(cr_mat)[i]
  ind1=abs(cr_mat[i,]) > 0.05 & abs(cr_mat[i,]) != 1
  cool=which(ind1, arr.ind = T)
  y_nam=rownames(as.data.frame(cool))
  print(x_nam) #Row name
  print(y_nam) #Column names with correlation >0.05
#}

  #Formula
form=as.formula(paste(x_nam,paste(y_nam,collapse = "+"),sep="~"))
glm(form,data=c1_data )

#GLMBoosting and Model Tuning, Depends on randomness
model1<-glmboost(form,data=c1_data,family = Gaussian(),
         center=TRUE,control = boost_control(mstop=200,nu=0.05,trace=TRUE))

#Induces randomness, can loop and take the nearest average integer
f<-cv(model1$`(weights)`,type="kfold",B=10)
cvm<-cvrisk(model1,folds=f)
opt_m<-mstop(cvm)

rmse<-function(model,col=i,data=c1_data){
  error=sqrt(mean((predict(model,data)-data[[col]])^2))
  return(error)
}

#Choosing the optimal model
model1[opt_m]
coef(model1,which="")
error=rmse(model1,col=i)
print(error)

