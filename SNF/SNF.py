#===================Importing Modules====================
from __future__ import division
import numpy
import pandas
#============Importing the Datasets======================
df1=pandas.read_csv("../bacterial_microbiome.csv")
df1.set_index("PatientID",inplace=True)
df2=pandas.read_csv("../fungal_microbiome.csv")
df2.set_index("PatientID",inplace=True)

print "Nan values present", df1.isnull().values.any()
print "Nan values present", df2.isnull().values.any()
print "length of df1", len(df1)
print "length of df2", len(df2)
#-------------Dropping NAN-------
df1.dropna(axis=0,how='any',inplace=True)
df2.dropna(axis=0,how='any',inplace=True)

print "length of df1", len(df1)
print "length of df2", len(df2)

#----------Making equal number of patients-------------
i1=set(df1.index)
i2=set(df2.index)
i=i1.intersection(i2)
print "length of intersection", len(i)

b_data=df1.loc[list(i),:]
f_data=df2.loc[list(i),:]
#============Bray Curtis Similarity Calculation ==========
from skbio.diversity.beta import pw_distances as dis
W1=(dis(b_data.values,ids=b_data.index).data - 1)*-1
W2=(dis(f_data.values,ids=f_data.index).data - 1)*-1
#============Calculation of P matrix=====================
def P(W):
    m=len(W)
    p=numpy.diag([0.5]*m)
    for i in range(m):
        for j in range(m):
            if j==i:
                pass
            else:
                p[i,j]=W[i,j]/(2*(numpy.sum(W[i,:])-W[i,i]))
    return p
#============Calculation of S matrix====================
'''
#S with k whole population

def S(W):
    m=len(W)
    s=numpy.zeros((m,m))
    for i in range(m):
        s[i,:]=W[i,:]/numpy.sum(W[i,:])
    return s
'''

# S with K nearest neighbours from snfpy
def S_k(W, K):
    m, n = W.shape
    IW1 = numpy.flip(W.argsort(axis=1), axis=1)
    newW = numpy.zeros(W.size)
    I1 = ((IW1[:, :K] * m) + numpy.vstack(numpy.arange(n))).flatten(order='F')
    newW[I1] = W.flatten(order='F')[I1]
    newW = newW.reshape(W.shape, order='F')
    newW = newW / newW.sum(axis=1)[:, numpy.newaxis]

    return newW
#============The merging step===========================
def SNF_merge(S_M,K,t):
    S_M=[W1,W2]
    P1=P(W1)
    P2=P(W2)
    S1=S_k(W1,229)
    S2=S_k(W2,229)
    iter=0
    while iter<t and numpy.sum(abs(P2-P1))>0.1:
        print "Iteration", iter
        print "distance",numpy.sum(abs(P2-P1))
        P1=S1.dot(P2).dot(S1.transpose())
        P1=P(P1)
        P2=S2.dot(P1).dot(S2.transpose())
        P2=P(P2)
        iter += 1        
    P_c=(P1+P2)/2
    return P_c
P_m=SNF_merge([W1,W2],229,10)

#================Clustering====================
from sklearn.decomposition import PCA
from sklearn.utils.validation import check_array

def get_n_clusters(arr, n_clusters=range(2, 10)):
    n_clusters = check_array(n_clusters, ensure_2d=False)
    eigenvalue = PCA().fit(arr).singular_values_[:-1]
    eigengap = numpy.abs(numpy.diff(eigenvalue))
    eigengap = eigengap * (1 - eigenvalue[:-1]) / (1 - eigenvalue[1:])
    n = eigengap[n_clusters - 1].argsort()[::-1]

    return n_clusters[n[0]], n_clusters[n[1]]

print "Optimal number of clusters using eigen gap", get_n_clusters(P_m)

from sklearn.cluster import SpectralClustering
sc = SpectralClustering(2, affinity='precomputed', n_init=100, assign_labels='discretize')
sc.fit(P_m)



def _silhouette_samples(arr, labels):
    from sklearn.preprocessing import LabelEncoder
    from sklearn.utils import check_X_y

    def check_number_of_labels(n_labels, n_samples):
        if not 1 < n_labels < n_samples:
            raise ValueError("Number of labels is %d. Valid values are 2 "
                             "to n_samples - 1 (inclusive)" % n_labels)

    arr, labels = check_X_y(arr, labels, accept_sparse=['csc', 'csr'])
    le = LabelEncoder()
    labels = le.fit_transform(labels)
    check_number_of_labels(len(le.classes_), arr.shape[0])

    unique_labels = le.classes_
    n_samples_per_label = numpy.bincount(labels, minlength=len(unique_labels))

    # For sample i, store the mean distance of the cluster to which
    # it belongs in intra_clust_dists[i]
    intra_clust_aff = numpy.zeros(arr.shape[0], dtype=arr.dtype)

    # For sample i, store the mean distance of the second closest
    # cluster in inter_clust_dists[i]
    inter_clust_aff = intra_clust_aff.copy()

    for curr_label in range(len(unique_labels)):

        # Find inter_clust_dist for all samples belonging to the same
        # label.
        mask = labels == curr_label
        current_distances = arr[mask]

        # Leave out current sample.
        n_samples_curr_lab = n_samples_per_label[curr_label] - 1
        if n_samples_curr_lab != 0:
            intra_clust_aff[mask] = numpy.sum(
                current_distances[:, mask], axis=1) / n_samples_curr_lab

        # Now iterate over all other labels, finding the mean
        # cluster distance that is closest to every sample.
        for other_label in range(len(unique_labels)):
            if other_label != curr_label:
                other_mask = labels == other_label
                other_distances = numpy.mean(
                    current_distances[:, other_mask], axis=1)
                inter_clust_aff[mask] = numpy.maximum(
                    inter_clust_aff[mask], other_distances)

    sil_samples = intra_clust_aff - inter_clust_aff
    sil_samples /= numpy.maximum(intra_clust_aff, inter_clust_aff)

    # score 0 for clusters of size 1, according to the paper
    sil_samples[n_samples_per_label.take(labels) == 1] = 0

    return sil_samples


def silhouette_score(arr, labels):
    return numpy.mean(_silhouette_samples(arr, labels))


print "Silhoutee score", silhouette_score(P_m,sc.labels_)
