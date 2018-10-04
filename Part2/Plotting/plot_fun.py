import pandas
import networkx as nx
import matplotlib.pyplot as plt
Adj=pandas.read_csv("Adjacency_matrix.csv")
p=pandas.read_csv("p_values.csv")
Adj.set_index("Unnamed: 0",inplace=True)
p.set_index("Unnamed: 0",inplace=True)

#--------Selecting p-values cutoff-------------
def p_select(Adj,p,p_cut_off,adj_cut_off):
    p.fillna(1,inplace=True)
    ind=p>p_cut_off
    Adj[ind]=0
    Adj[abs(Adj)<adj_cut_off]=0
    return(Adj)
#-------------------------------------------
p.fillna(1,inplace=True)
p_cut_off=0.001
ind=p>p_cut_off
Adj[ind]=0

#Adjacency matrix cut-off
Adj[abs(Adj)<750]=0

#---------Histogram of edge weights--------
plt.hist(Adj.values.flatten(),bins=50)
plt.show()
#-----------------------------------------------

labels={}
col=Adj.columns.values
for i in range(len(col)):
    labels[i]=str(col[i])

G=nx.DiGraph(Adj.values)

to_rem=list(nx.isolates(G))
G.remove_nodes_from(to_rem)
map(labels.pop,to_rem)
nx.draw(G,labels=labels)
plt.show()
