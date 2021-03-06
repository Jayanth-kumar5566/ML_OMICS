library(dunn.test)
#Combining the data
data1=read.csv("Documents/MS/ML_OMICS/SNF_standard/Cluster_eval/DM_CM.csv")
data2=read.csv("Documents/MS/ML_OMICS/SNF_standard/labels.csv")
#Add patient ID in the labels file manually(from luminex.csv)
data=merge(data1,data2,by="Patient_ID")

rownames(data)<-data[["Patient_ID"]]
data=data[-1]


  
#data=na.omit(data) #Add NA in the dataset for missing values(1 missing present)

#Splitting the data based clusters
X=split(data,data$x)

pop1=X$'1'
pop2=X$'2'
pop3=X$`3`
pop4=X$`4`
pop5=X$`5`

nam=colnames(data)

sink("Documents/MS/ML_OMICS/SNF_standard/4_Clust_SNF/post_hoc.txt") 
#--------------For 4 Clusters Outcome Comparision-----------------------

#Numerical
for (i in nam[1:8])
{
p1=na.omit(pop1[[i]])
p2=na.omit(pop2[[i]])
p3=na.omit(pop3[[i]])
p4=na.omit(pop4[[i]])
boxplot(p1,p2,p3,p4,xlab=i)

#p1=na.omit(p1[!p1 %in% boxplot.stats(p1)$out])
#p2=na.omit(p2[!p2 %in% boxplot.stats(p2)$out])
#p3=na.omit(p3[!p3 %in% boxplot.stats(p3)$out])

t=kruskal.test(x=list(p1,p2,p3,p4))
#t=wilcox.test(p1,p2,alternative = "two.sided",paired=FALSE)
if (t$p.value <= 0.05){
  print(i)
  d=dunn.test(x=list(p1,p2,p3,p4),method = 'bh',alpha=0.05,kw=FALSE,list = TRUE,table = FALSE,altp = TRUE)
  # png(paste("./Documents/MS/ML_OMICS/SNF_4/Cluster_eval/",i,".png",sep=''))
  #boxplot(p1,p2,p3,xlab=i)
  # dev.off()
  #print(c(median(p4,na.rm = TRUE),"+/-",var(p4,na.rm=TRUE)))
}
}

for (i in nam[9:16]){
  le=levels(factor(data[,i]))
  pop1[i]<-lapply(pop1[i],factor,levels=le)
  pop2[i]<-lapply(pop2[i],factor,levels=le)
  pop3[i]<-lapply(pop3[i],factor,levels=le)
  pop4[i]<-lapply(pop4[i],factor,levels=le)
  pop5[i]<-lapply(pop5[i],factor,levels=le)
}

#Categorical
for (i in nam[9:16])
{
p1=pop1[[i]]
p2=pop2[[i]]
p3=pop3[[i]]
p4=pop4[[i]]

p1=na.omit(p1)
p2=na.omit(p2)
p3=na.omit(p3)
p4=na.omit(p4)

x=rbind(summary(p1),summary(p2),summary(p3),summary(p4))
#print(i)
#print(x)
#t=fisher.test(x,simulate.p.value = TRUE)
t=chisq.test(x,simulate.p.value = TRUE)
#print(t$p.value)
if (t$p.value <= 0.05){
  print(i)
  d=dunn.test(x=list(p1,p2,p3,p4),method = 'bh',alpha=0.05,kw=FALSE,list = TRUE,table = FALSE,altp = TRUE)
  # for (i in 1:4){
  #   print(x[i,]/sum(x[i,]))
  # }
  # png(paste("./Documents/MS/ML_OMICS/SNF_4/Cluster_eval/",i,".png",sep=''))
  par(mfrow=c(2,2))
  li=max(x)
  plot(p1,xlab=i,ylim=c(1,li))
  plot(p2,xlab=i,ylim=c(1,li))
  plot(p3,xlab=i,ylim=c(1,li))
  plot(p4,xlab = i,ylim = c(1,li))
  # dev.off()
  }
}
sink()
#------------------Median or mean values for each population----------------

sink("Documents/MS/ML_OMICS/SNF_standard/4_Clust_SNF/numerical.txt")
#Numerical
for (i in nam[1:8])
{
  p1=pop1[[i]]
  p2=pop2[[i]]
  p3=pop3[[i]]
  p4=pop4[[i]]
 
  t=kruskal.test(list(p1,p2,p3,p4))
  if (t$p.value <= 0.05){
    print(i)
    count=1
    for (j in list(p1,p2,p3,p4)){
      print(paste("C",count,collapse = ""))
      if (shapiro.test(j)$p.value<0.05){
        print("Median +/- variance")
        print(paste(median(j,na.rm = TRUE),"+/-",var(j,na.rm=TRUE),collapse=""))
      }
      else{
        print("Mean +/- std")
      print(paste(mean(j,na.rm = TRUE),"+/-",sqrt(var(j,na.rm=TRUE)),collapse=""))
      }
      count = count +1
    }
  }
}

sink()
#---------------------Categorical--------------------------
sink("Documents/MS/ML_OMICS/SNF_standard/4_Clust_SNF/categ_toprocess.txt")

for (i in nam[9:16]){
  le=levels(factor(data[,i]))
  pop1[i]<-lapply(pop1[i],factor,levels=le)
  pop2[i]<-lapply(pop2[i],factor,levels=le)
  pop3[i]<-lapply(pop3[i],factor,levels=le)
  pop4[i]<-lapply(pop4[i],factor,levels=le)
}

#Categorical
for (i in nam[9:16])
{
  p1=pop1[[i]]
  p2=pop2[[i]]
  p3=pop3[[i]]
  p4=pop4[[i]]


  p1=na.omit(p1)
  p2=na.omit(p2)
  p3=na.omit(p3)
  p4=na.omit(p4)

  x=rbind(summary(p1),summary(p2),summary(p3),summary(p4))
  t=chisq.test(x,simulate.p.value = TRUE)
  if (t$p.value <= 0.05){
    print(i)
    print((x/rowSums(x))*100)
    }
  }
sink()
