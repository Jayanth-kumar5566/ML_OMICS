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

#Adjacency matrix and p-value matrix creation
m=matrix(0,ncol=dim(cr_mat)[1],nrow=dim(cr_mat)[2])
p_m=matrix(2,ncol=dim(cr_mat)[1],nrow=dim(cr_mat)[2])
m<-data.frame(m,row.names = colnames(c1_data))
p_m<-data.frame(p_m,row.names = colnames(c1_data))
colnames(m)<-colnames(c1_data)
colnames(p_m)<-colnames(c1_data)

r2<-function(model,col=i,data=c1_data){
  sse=sum((predict(model,data)-data[[col]])^2)
  tss=sum((data[[col]]-mean(data[[col]]))^2)
  error=1-(sse/tss)
  return(error)
}

for (i in 1:539){
  x_nam=rownames(cr_mat)[i]
  ind1=abs(cr_mat[i,]) > 0.05 & abs(cr_mat[i,]) != 1
  cool=which(ind1, arr.ind = T)
  y_nam=rownames(as.data.frame(cool)) #Column names with correlation >0.05
  if(identical(y_nam,character(0)) == TRUE){next}
  print(x_nam) #Row name
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
  if(opt_m==0){opt_m=1}
  #Choosing the optimal model
  model1[opt_m]
  error=r2(model1,col=i)
  if (error<0.5){next}
  wghts<-coef(model1,which="")
  x<-t(as.data.frame(wghts[-1]))
  row.names(x)<-x_nam
    #Appending the coefficient matrix to adjacency matrix
    for(cl in colnames(x)){
      m[x_nam,cl]<-x[x_nam,cl]
    }
  #Bootstrap distribution
  #------------Using boot-----------------  
  library(boot)
  boot.stat<-function(data,indices,m_stop,form,x_nam){
    data<-data[indices,]
    mod<-glmboost(form,data=data,family = Gaussian(),
                  center=FALSE,control = boost_control(mstop=m_stop,nu=0.05,trace=TRUE))
    wghts<-coef(mod,which="")
    x<-t(as.data.frame(wghts[-1]))
    row.names(x)<-x_nam
    return(x)
  }
  
  model.boot<-boot(c1_data,boot.stat,100,m_stop=opt_m,form=form,x_nam=x_nam)
  #-----------Permutation with renormalization-------------
  #copy of the data
  c1_data_p=c1_data
  out_comb<-x
  #permutation
  counter=0
  while (counter<100){
    c1_data_p[[x_nam]]<-sample(c1_data_p[[x_nam]])
    #renormalization
    c1_data_p=(c1_data_p/rowSums(c1_data_p))*200
    out<-boot.stat(c1_data_p,indices = 1:dim(c1_data)[1],m_stop=opt_m,form=form,x_nam=x_nam)
    out_comb=rbind(out_comb,out)
    counter = counter + 1
  }
  out_comb<-out_comb[-1,]
  #Comparing two distributions
  for (i in 1:dim(out_comb)[2]){
    p=wilcox.test(model.boot$t[,i],out_comb[,i],alternative = "two.sided",paired = FALSE)$p.value
    p_m[x_nam,colnames(out_comb)[i]]<-p
  }
}

#To do correction of multiple comparision later

#model.boot$t[,1] #hist values for 100 iterations of first parameter

#out_comb[,1] #values of the first parameter
#--------------------------------------------------------
