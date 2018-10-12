import pandas
import networkx as nx
import matplotlib.pyplot as plt
Adj=pandas.read_csv("Adjacency_matrix.csv")
p=pandas.read_csv("p_values.csv")
ab=pandas.read_csv("c2_abund.csv")
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

#---------Histogram of edge weights--------
Ad_f=Adj.values.flatten()
plt.hist(Ad_f)
plt.show()

plt.boxplot(Ad_f)
plt.show()
#Adj[abs(Adj)<750]=0
#-----------------------------------------------

Ad_f.sort()
#Adjacency matrix cut-off
'''
#Top 20% and Below 20%
bel=int(round(len(Ad_f)*0.0001)) 
top=int(round(len(Ad_f)*0.9999))
#Top 20 and below 20 based on the value 
bel=int(round(max(Ad_f)*0.20))
top=int(round(max(Ad_f)*0.20))

#Adjacency matrix cut-off
Adj[(Ad_f[bel]<Adj) & (Adj<Ad_f[top])]=0
'''
top=float(raw_input("Top cutoff"))
bel=float(raw_input("Below cutoff"))
Adj[(bel<Adj) & (Adj<top)]=0

#----------------------------------------------

labels={}
col=Adj.columns.values
for i in range(len(col)):
    labels[i]=str(col[i])

G=nx.DiGraph(Adj.values)
to_rem=list(nx.isolates(G))
G.remove_nodes_from(to_rem)
map(labels.pop,to_rem)

ns=[ab.ix[i,'x'] for i in G.nodes()]

edges,weights = zip(*nx.get_edge_attributes(G,'weight').items())
pos=nx.random_layout(G)
nx.draw(G,labels=labels,font_size=10,node_size=ns,edgelist=edges,edge_color=weights,edge_cmap=plt.cm.seismic,node_color='r',arrows=False,width=3)
plt.show()
