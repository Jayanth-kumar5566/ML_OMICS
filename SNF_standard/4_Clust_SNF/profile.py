import pandas
import matplotlib.pyplot as plt
import numpy
from scipy.stats.mstats import kruskalwallis

df=pandas.read_csv("fungi.csv")
df.dropna(axis=0,how="any",inplace=True)
df.set_index("PatientID",inplace=True)
pop1=df[df.x==1]
pop2=df[df.x==2]
pop3=df[df.x==3]
pop4=df[df.x==4]
'''
pop3=df[df.x==3]
pop4=df[df.x==4]
pop5=df[df.x==5]
'''

mat=pop4.values[:,:-1]
x=pop4.columns[:-1]

'''
from mpl_toolkits.axes_grid1 import make_axes_locatable
ax = plt.subplot(111)
ax.set_xticks(range(mat.shape[1]))
ax.set_xticklabels(x,rotation='vertical',size=10)
im = ax.imshow(mat,cmap="Greys",aspect='auto')
plt.title("Fungi Cluster 1")
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.show()


'''

#Plot average
y=numpy.median(mat,axis=0)
y=y.reshape(1,y.shape[0])
from mpl_toolkits.axes_grid1 import make_axes_locatable
ax = plt.subplot(111)
ax.set_xticks(range(y.shape[1]))
ax.set_xticklabels(x,rotation='vertical',size=10)
im = ax.imshow(y,cmap="Reds",aspect='auto')
plt.title("Bacterial Cluster 1 average")
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
#plt.show()

d=dict()
for i in range(len(x)):
    d[str(x[i])]=y[0][i]
'''    
import operator
d=sorted(d.items(),key=operator.itemgetter(1))

file=open("/tmp/x.txt",'w')
for i in d[::-1]:
    file.write(str(i[0])+": "+str(i[1])+"\n")

file.close()
'''
#------------------Maan Whitney U test---------------------
file=open("/tmp/y.csv",'w')
for i in x:
    if sum(pop1[i].values)+sum(pop3[i].values)==sum(pop2[i].values)+sum(pop4[i].values)==0:
        pass
    else:
        (stat,p_value)=kruskalwallis(pop1[i].values,pop2[i].values,pop3[i].values,pop4[i].values)
        #Testing using 2 sided maan whitney for equal means
        if p_value <=0.05:
            print i
            file.write(str(i)+","+str(d[i])+","+str(p_value)+"\n")
file.close()
            
