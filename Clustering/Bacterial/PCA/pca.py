import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas

from sklearn import decomposition
from sklearn import datasets

X=pandas.read_csv("../bacterial_microbiome.csv")
X.dropna(axis=0,how='any',inplace=True)


X=X.replace(to_replace="Mild",value=0)
X=X.replace(to_replace="Moderate",value=1)
X=X.replace(to_replace="Severe",value=2)

#np.random.seed(5)
X=X.as_matrix()
#Y=Y.as_matrix()

y1=X[:,61].astype('int') #BSI Score
y2=X[:,60].astype('int') #Mild severe and Moderate
X=X[:,1:60].astype('int') #Data


fig = plt.figure(1, figsize=(4, 3))
plt.clf()
ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

plt.cla()
pca = decomposition.SparsePCA(n_components=3)
pca.fit(X)
X = pca.transform(X)

ax.scatter(X[:, 0], X[:, 1],X[:, 2],c=y2,s=50)

ax.w_xaxis.set_ticklabels([])
ax.w_yaxis.set_ticklabels([])
ax.w_zaxis.set_ticklabels([])

plt.show()


'''
#==========2D PCA====================
pca = decomposition.SparsePCA(n_components=3)
pca.fit(X)
X = pca.transform(X)

plt.scatter(X[:, 0], X[:, 1],c=y2,s=50)

plt.show()

'''
