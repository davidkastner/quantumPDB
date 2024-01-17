import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def format_plot() -> None:
    """
    General plotting parameters for the Kulik Lab.
    """
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

# Load the CSV files into pandas DataFrames
no_smoothing = pd.read_csv('no_smoothing.csv')
box_plot = pd.read_csv('quartile.csv')
dbscan = pd.read_csv('dbscan.csv')
dummy = pd.read_csv('dummy.csv')

# Label the data
no_smoothing['Method'] = 'none'
box_plot['Method'] = 'quartile'
dbscan['Method'] = 'dbscan'
dummy['Method'] = 'dummy'

# Concatenate all data into a single DataFrame
all_data = pd.concat([no_smoothing, box_plot, dbscan, dummy])

# Select only the relevant columns (assuming the values are in the third column)
all_data = all_data.iloc[:, [2, -1]]
all_data.columns = ['Value', 'Method']

colors = ["#219ebc", "#023047", "#ffb703", "#fb8500"]

# Create the strip plot
format_plot()
plt.figure(figsize=(6, 4))
plt.ylim(0, 5000)
sns.stripplot(x='Method', y='Value', data=all_data, jitter=True, palette=colors, alpha=.5)

# Calculate medians for each method and plot horizontal lines
medians = all_data.groupby('Method')['Value'].median()
methods = all_data['Method'].unique()

for method in methods:
    median_val = medians[method]
    plt.hlines(median_val,xmin=methods.tolist().index(method)-0.12, xmax=methods.tolist().index(method)+0.12, color='black', lw=4, zorder=5, alpha=1)

plt.xlabel('smoothing method', weight="bold")
plt.ylabel('atom count', weight="bold")
plt.savefig("strip_plot.png", bbox_inches="tight", format="png", dpi=600)