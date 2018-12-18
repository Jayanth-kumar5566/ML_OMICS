import matplotlib.pyplot as plt
import seaborn as sns
import pandas
M=pandas.read_csv("merged_matrix.csv")
lab=pandas.read_csv("labels.csv")
del M["Unnamed: 0"]
M.set_index(lab["Patient_ID"],inplace=True)
M.columns=lab["Patient_ID"]
lab.sort_values('x',inplace=True)
M=M.reindex(lab['Patient_ID'],axis='index')
M=M.reindex(lab['Patient_ID'],axis='columns')
ax=sns.heatmap(M,xticklabels=False,yticklabels=False,cmap=plt.get_cmap("Blues"))
plt.show()
