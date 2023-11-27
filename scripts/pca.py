import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import numpy as np

def format_plot() -> None:
    """General plotting parameters for the Kulik Lab."""

    font = {"family": "sans-serif", "weight": "bold", "size": 10}
    plt.rc("font", **font)
    plt.rcParams["xtick.major.pad"] = 5
    plt.rcParams["ytick.major.pad"] = 5
    plt.rcParams["axes.linewidth"] = 2
    plt.rcParams["xtick.major.size"] = 7
    plt.rcParams["xtick.major.width"] = 2
    plt.rcParams["ytick.major.size"] = 7
    plt.rcParams["ytick.major.width"] = 2
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["svg.fonttype"] = "none"

format_plot()

# Load your dataset
file_path = 'charge_results.csv'
data = pd.read_csv(file_path)

# Selecting the features for PCA
features = [
    '1_neg_charge_count', '1_pos_charge_count', '1_neg_charge_density', '1_pos_charge_density',
    '2_neg_charge_count', '2_pos_charge_count', '2_neg_charge_density', '2_pos_charge_density',
    '3_neg_charge_count', '3_pos_charge_count', '3_neg_charge_density', '3_pos_charge_density'
]
x = data.loc[:, features].values
pdb_ids = data['pdb_id'].values  # Assuming 'pdb_id' is the column name

# Standardizing the features
x = StandardScaler().fit_transform(x)

# Performing PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])

# Print the composition of principal components with feature names
print("PCA Component 1:")
for weight, feature in sorted(zip(pca.components_[0], features), reverse=True):
    print(f"{feature}: {weight}")

print("\nPCA Component 2:")
for weight, feature in sorted(zip(pca.components_[1], features), reverse=True):
    print(f"{feature}: {weight}")

# Plotting the results
plt.figure()
plt.axhline(y=0, color='k', linestyle='-', linewidth=2, alpha=.2, zorder=1)
plt.axvline(x=0, color='k', linestyle='-', linewidth=2, alpha=.2, zorder=1)
plt.xlabel('Principal Component 1', fontsize=12, fontweight='bold')
plt.ylabel('Principal Component 2', fontsize=12, fontweight='bold')
plt.scatter(principalDf['PC1'], principalDf['PC2'], alpha=0.25, color="blue")

# Labeling outliers
distances = np.sqrt(principalDf['PC1']**2 + principalDf['PC2']**2)
outliers = distances.nlargest(8).index  # Label the 5 farthest points; adjust number as needed
for outlier in outliers:
    plt.annotate(pdb_ids[outlier], (principalDf.iloc[outlier]['PC1'], principalDf.iloc[outlier]['PC2']),
                 fontweight='bold')

# Saving the plot as a PNG file
plt.savefig('pca_result_with_outliers.png')
