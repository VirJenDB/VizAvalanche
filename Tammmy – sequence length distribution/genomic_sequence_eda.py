import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import seaborn as sns

DATASET_PATH = '/home/srichter/Schreibtisch/plot_data/collected_data.csv'

print(f"Loading dataset from {DATASET_PATH}...")
try:
    df = pd.read_csv(DATASET_PATH)
    print(f"Dataset loaded successfully. Shape: {df.shape}")
    print("\nFirst 5 rows:")
    print(df.head())
    print("\nDataset info:")
    df.info()
except FileNotFoundError:
    print(f"Error: The file '{DATASET_PATH}' was not found. Please ensure it was unzipped correctly and the path is accurate.")
except Exception as e:
    print(f"An error occurred while loading the dataset: {e}")

test_size_ratio = 0.2

df_train, df_test = train_test_split(df, test_size=test_size_ratio, random_state=42)

print(f"Original dataset shape: {df.shape}")
print(f"Training set shape: {df_train.shape}")
print(f"Test set shape: {df_test.shape}")

print("\nFirst 5 rows of the training set:")
print(df_train.head())

print("\nFirst 5 rows of the test set:")
print(df_test.head())

print("Generating sequence length distribution visualization")
plt.figure(figsize=(12, 7))
sns.histplot(df['sequence_length'], bins=50, kde=True, log_scale=True)
plt.title('Sequence Length Distribution (Log Scale X & Y)', fontsize=16)
plt.xlabel('Sequence Length (Log Scale)', fontsize=12)
plt.ylabel('Count (Log Scale)', fontsize=12)
plt.yscale('log') # Add log scale to the Y-axis
plt.grid(axis='both', alpha=0.75) # Grid for both axes
plt.tight_layout()
plt.show()
