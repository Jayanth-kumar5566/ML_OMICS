from __future__ import division
import pandas
import numpy
df=pandas.read_csv("luminex.csv")
df.set_index("PatientID",inplace=True)
def normalize(x):
    t=numpy.sum(x)
    for i in range(len(x)):
        x[i]=(x[i]/t)*100
    return x
df.apply(normalize,axis=1)
df.to_csv("luminex.csv")
