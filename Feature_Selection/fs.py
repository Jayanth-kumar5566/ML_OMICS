import pandas
df=pandas.read_csv("../bacterial_microbiome.csv")

#====Finding the Maximum of the columns=========

df_max=df.max().to_dict()
cutoff=1 #Set the cutoff here
cols=[] #contains all the filtered columns

for i in df_max:
    if df_max[str(i)]>=1:
        cols.append(str(i))
    else:
        pass

df_filtered=df[cols]

#==========Save it as a csv file===========
df_filtered=df_filtered.set_index("PatientID")
df_filtered.to_csv("./Selected/bacterial_microbiome.csv")

