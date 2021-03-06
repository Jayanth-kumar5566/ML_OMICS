import pandas
import matplotlib.pyplot as plt
import numpy
from scipy.stats import mannwhitneyu

df=pandas.read_csv("fungi.csv")
df.dropna(axis=0,how="any",inplace=True)
df.set_index("PatientID",inplace=True)
pop1=df[df.x==1]
pop2=df[df.x==2]

'''
pop3=df[df.x==3]
pop4=df[df.x==4]
pop5=df[df.x==5]
'''

df1=pop2

df1=(df1.T[df1.sum(axis=0)!=0]).T
mat=df1.values[:,:-1]
x=df1.columns[:-1]

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
y_v=numpy.var(mat,axis=0)
'''
y=y.reshape(1,y.shape[0])
y_v=y_v.reshape(1,y_v.shape[0])
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
'''
d={"species":x,"median":y,"variance":y_v}
d=pandas.DataFrame(data=d)
d.sort_values("variance",ascending=True,axis=0,inplace=True)
d.to_csv("/tmp/fungi_c2.csv")

'''
#------------------Maan Whitney U test---------------------
file=open("/tmp/y.csv",'w')
for i in x:
    if sum(pop1[i].values)==sum(pop2[i].values)==0:
        pass
    else:
        (stat,p_value)=mannwhitneyu(pop1[i],pop2[i])
        #Testing using 2 sided maan whitney for equal means
        if p_value*2 <=0.05:
            print i
            file.write(str(i)+","+str(d[i])+","+str(p_value)+"\n")
file.close()
            
'''
