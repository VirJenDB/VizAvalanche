import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import matplotlib.pyplot as plt
import seaborn as sns



df = pd.read_csv("is_cluster_representative.csv")
# print(df.columns)


df['collection_date'] = pd.to_datetime(df['collection_date'], errors='coerce')
submissions = df.groupby(df['collection_date'].dt.to_period("Y")).size()

submissions.plot(kind='bar')
plt.xlabel("Year")
plt.ylabel("Number of submissions")
plt.title("Submissions over time")
plt.savefig("Submissions over time.png")
plt.show()







df2 = pd.read_csv("sample_5k.csv")





df2['sequence_gc_content'].hist(bins=50)
plt.xlabel("GC content (%)")
plt.ylabel("Number of sequences")
plt.title("Distribution of GC content")
plt.savefig("Distribution of GC content.png")
plt.show()







plt.boxplot(df2['sequence_gc_content'].dropna())
plt.ylabel("GC content (%)")
plt.title("GC content spread")
plt.savefig("GC content spread.png")
plt.show()






sns.boxplot(x=df['molecule_type'], y=df2['sequence_gc_content'])
plt.title("GC content by molecule type")
plt.savefig("GC content by molecule type.png")
plt.show()




import pandas as pd
import matplotlib.pyplot as plt

# قراءة الملفين
df_cluster = pd.read_csv("is_cluster_representative.csv", usecols=['sequence_length','collection_date'])
df_sample = pd.read_csv("sample_5k.csv", usecols=['sequence_length','sequence_gc_content'])

# دمج باستخدام sequence_length
merged = pd.merge(df_cluster, df_sample, on='sequence_length', how='inner')

# تحويل التاريخ
merged['collection_date'] = pd.to_datetime(merged['collection_date'], errors='coerce')

# حساب متوسط GC لكل سنة
gc_trend = merged.groupby(merged['collection_date'].dt.to_period("Y"))['sequence_gc_content'].mean()

# رسم الخط
gc_trend.plot(kind='line', marker='o')
plt.xlabel("Year")
plt.ylabel("Average GC content (%)")
plt.title("GC content trend over time")
plt.savefig("GC content trend over time.png")
plt.show()
