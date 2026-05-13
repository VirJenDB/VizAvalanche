import pandas as pd
import plotly.express as px

# 1. LOAD YOUR DATA
# Replace 'your_data.csv' with the actual path to the file you downloaded
df = pd.read_csv('/Users/tamanasisodiya/Downloads/collected_data.csv')

# --- IDEA 1: Sequence Length Distribution ---
def plot_sequence_length(data):
    # 'length' is a common field name, check your CSV columns if it's different
    fig = px.histogram(data, x="length", 
                       title="Sequence Length Distribution",
                       labels={'length': 'Sequence Length (bp)'},
                       nbins=100,
                       log_y=True) # Adding log scale as requested in your ideas!
    fig.show()

# --- IDEA 2: GC Content Bins ---
def plot_gc_content(data):
    # 'gc_content' is usually the field name
    fig = px.histogram(data, x="gc_content", 
                       title="GC Content Distribution",
                       nbins=50,
                       color_discrete_sequence=['indianred'])
    fig.show()

# Run the functions
plot_sequence_length(df)
plot_gc_content(df)