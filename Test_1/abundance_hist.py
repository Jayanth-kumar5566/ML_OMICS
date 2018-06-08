import matplotlib.pyplot as plt
import pandas

df=pandas.read_csv("./bacterial_microbiome.csv")
df.dropna(axis=0,how='any',inplace=True)
df.set_index("PatientID",inplace=True)
patno=df.index.tolist()

patients=df[df.Status=="Bronchiectasis"]
controls=df[df.Status=="Healthy"]

del patients["Status"]
del controls["Status"]


x_ticks=patients.columns.tolist()

plt.xticks(range(59), x_ticks, rotation=90)

for i in patno:
    if i in patients.index.tolist():
       plt.scatter(range(59),patients.loc[str(i)],color='red')
    elif i in controls.index.tolist():
        plt.scatter(range(59),controls.loc[str(i)],color='blue')


plt.show()


#=============Plotting the mean of them====================

plt.xticks(range(59), x_ticks, rotation=90,size=8)
plt.scatter(range(59),patients.mean(),color='red')
plt.scatter(range(59),controls.mean(),color='blue')
plt.show()

    
