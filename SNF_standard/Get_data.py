#===================Importing Modules====================
from __future__ import division
import numpy
import pandas
#============Importing the Datasets======================
df1=pandas.read_csv("../bacterial_microbiome.csv")
df1.set_index("PatientID",inplace=True)
df2=pandas.read_csv("../fungal_microbiome.csv")
df2.set_index("PatientID",inplace=True)
df3=pandas.read_csv("../Luminex.csv")
df3.set_index("PatientID",inplace=True)

print "Nan values present", df1.isnull().values.any()
print "Nan values present", df2.isnull().values.any()
#print "Nan values present", df3.isnull().values.any()
print "length of df1", len(df1)
print "length of df2", len(df2)
#print "length of df3", len(df3)
#-------------Dropping NAN-------
df1.dropna(axis=0,how='any',inplace=True)
df2.dropna(axis=0,how='any',inplace=True)
#df3.dropna(axis=0,how='any',inplace=True)

print "length of df1", len(df1)
print "length of df2", len(df2)
#print "length of df2", len(df3)
#----------Making equal number of patients-------------
i1=set(df1.index)
i2=set(df2.index)
#i3=set(df3.index)
i=i1.intersection(i2)
#i=i.intersection(i3)
print "length of intersection", len(i)

b_data=df1.loc[list(i),:]
f_data=df2.loc[list(i),:]
l_data=df3.loc[list(i),:]

#l_data.dropna(axis=0,how="any",inplace=True)

b_data.to_csv("bacteria.csv")
f_data.to_csv("fungi.csv")
l_data.to_csv("luminex.csv")
