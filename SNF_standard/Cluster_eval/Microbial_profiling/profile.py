import pandas
import matplotlib.pyplot as plt
import numpy
df=pandas.read_csv("luminex.csv")
df.set_index("PatientID",inplace=True)
pop1=df[df.x==1]
pop2=df[df.x==2]
pop3=df[df.x==3]
pop4=df[df.x==4]
pop5=df[df.x==5]


mat=pop5.values[:,:-1]
x=pop5.columns[:-1]

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
import operator
d=sorted(d.items(),key=operator.itemgetter(1))

file=open("/tmp/x.txt",'w')
for i in d[::-1]:
    file.write(str(i[0])+": "+str(i[1])+"\n")

file.close()

