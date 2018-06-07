data=read.csv("./../bacterial_microbiome.csv")
labels=read.csv("./../../../Disease_Measures.csv")
labels=labels$Disease_severity
data[61]=labels
data=na.omit(data)
X=data[,2:60]
y=data[61]
x=as.matrix(X,type='numeric')
y=as.matrix(y,type='numeric')
require(skmeans)
hparty <- skmeans(x, 3,m=1, control = list(verbose = TRUE))
table(y,hparty$cluster)
