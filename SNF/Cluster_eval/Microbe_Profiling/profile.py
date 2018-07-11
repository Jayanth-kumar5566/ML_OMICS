import pandas
import matplotlib.pyplot as plt
df=pandas.read_csv("bacteria.csv")
df.set_index("PatientID",inplace=True)
pop1=df[df.x==1]
pop2=df[df.x==2]

mat=pop1.values[:,:-1]
x=pop1.columns[:-1]
'''
from mpl_toolkits.axes_grid1 import make_axes_locatable
ax = plt.subplot(111)
ax.set_xticks(range(mat.shape[1]))
ax.set_xticklabels(x,rotation='vertical',size=5)
im = ax.imshow(mat,cmap="Greys",aspect='auto')
plt.title("Fungi Cluster 1")
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.show()
'''


'''
x=range(mat.shape[1])
for i in range(mat.shape[0]):
    plt.plot(x,mat[i,:],'o')

plt.show()
'''

#Plot average
y=mat.mean(axis=0)
y=y.reshape(1,y.shape[0])
from mpl_toolkits.axes_grid1 import make_axes_locatable
ax = plt.subplot(111)
ax.set_xticks(range(y.shape[1]))
ax.set_xticklabels(x,rotation='vertical',size=5)
im = ax.imshow(y,cmap="Reds",aspect='auto')
plt.title("Bacterial Cluster 1 average")
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.show()
