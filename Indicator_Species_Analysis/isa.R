library("indicspecies")
data=read.csv("Documents/MS/ML_OMICS/Part2/data_combined.csv")
Xdata=data[,-541]
ydata=data[,541]
row.names(Xdata)<-Xdata$PatientID
Xdata$PatientID<-NULL
indval=multipatt(Xdata,ydata,control = how(nperm=999))
summary(indval)
summary(indval,indvalcomp = TRUE,alpha = 0.05)
indval$sign
