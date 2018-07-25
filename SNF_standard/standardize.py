from __future__ import division
import pandas
import numpy
from scipy.stats import zscore
df=pandas.read_csv("luminex.csv")
df1=pandas.read_csv("bacteria.csv")
df2=pandas.read_csv("fungi.csv")
df.set_index("PatientID",inplace=True)
df1.set_index("PatientID",inplace=True)
df2.set_index("PatientID",inplace=True)

df=df.loc[:, (df != 0).any(axis=0)]
df1=df1.loc[:, (df1 != 0).any(axis=0)]
df2=df2.loc[:, (df2 != 0).any(axis=0)]

def normalize(x):
    t=numpy.sum(x)
    for i in range(len(x)):
        x[i]=(x[i]/t)*100
    return x
'''
df=df.apply(zscore,axis=0)
df1=df1.apply(zscore,axis=0)
df2=df2.apply(zscore,axis=0)


#Translate the values such that all the values are positive
mini=numpy.min([numpy.min(df.values),numpy.min(df1.values),numpy.min(df2.values)])
df=df+abs(mini)
df1=df1+abs(mini)
df2=df2+abs(mini)
print "Setting the mean to ",mini
'''

df=df.apply(normalize,axis=1)
df1=df1.apply(normalize,axis=1)
df2=df2.apply(normalize,axis=1)

df.to_csv("luminex.csv")
df1.to_csv("bacteria.csv")
df2.to_csv("fungi.csv")
