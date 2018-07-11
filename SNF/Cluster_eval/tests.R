#Combining the data
data1=read.csv("Documents/MS/ML_OMICS/SNF/Cluster_eval/DM_CM.csv")
data2=read.csv("Documents/MS/ML_OMICS/SNF/Cluster_eval/labels.csv")

data=merge(data1,data2,by="Patient_ID")

rownames(data)<-data[["Patient_ID"]]
data=data[-1]


  
#data=na.omit(data) #Add NA in the dataset for missing values(1 missing present)

#Splitting the data based clusters
X=split(data,data$x)

pop1=X$'1'
pop2=X$'2'

nam=colnames(data)


#--------------For 3 Clusters Outcome Comparision-----------------------

#Numerical
for (i in nam[1:8])
{
p1=pop1[[i]]
p2=pop2[[i]]

#t=kruskal.test(list(p1,p2))
t=wilcox.test(p1,p2,alternative = "two.sided",paired=FALSE)
#print(t$p.value)
if (t$p.value <= 0.05){
  print(i)
  png(paste("./Documents/MS/ML_OMICS/SNF/Cluster_eval/",i,".png",sep=''))
  boxplot(p1,p2,xlab=i)
  dev.off()
}
}

for (i in nam[9:16]){
  pop1[i]<-lapply(pop1[i],factor)
  pop2[i]<-lapply(pop2[i],factor)
}

#Categorical
for (i in nam[9:16])
{
#i=10
p1=pop1[[i]]
p2=pop2[[i]]
p1=na.omit(p1)
p2=na.omit(p2)
x=rbind(summary(p1),summary(p2))
#print(i)
#print(x)
t=fisher.test(x,simulate.p.value = TRUE)
#t=chisq.test(x)
if (t$p.value <= 0.05){print(i)
  png(paste("./Documents/MS/ML_OMICS/SNF/Cluster_eval/",i,".png",sep=''))
  par(mfrow=c(1,2))
  li=max(x)
  plot(p1,xlab=i,ylim=c(1,li))
  plot(p2,xlab=i,ylim=c(1,li))
  dev.off()}
}

