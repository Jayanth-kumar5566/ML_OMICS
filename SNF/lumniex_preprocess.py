from __future__ import division
import pandas
import numpy
from scipy.stats import zscore
df=pandas.read_csv("luminex.csv")
df.set_index("PatientID",inplace=True)
def normalize(x):
    t=numpy.sum(x)
    for i in range(len(x)):
        x[i]=(x[i]/t)*100
    return x
df=df.apply(zscore,axis=0)
#Translate the values such that all the values are positive
mini=numpy.min(df.values)
df=df+abs(mini)

df=df.apply(normalize,axis=1)
df.to_csv("luminex.csv")
