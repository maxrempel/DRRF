import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.ticker import MaxNLocator

### Load all datasets
df1 = pd.read_csv("dataset1_all_chromosomes_exbor1.tsv", sep="\t")
df2 = pd.read_csv("dataset1_all_chromosomes_exbor2.tsv", sep="\t")
df3 = pd.read_csv("dataset2_all_chromosomes_exbor1.tsv", sep="\t")
df4 = pd.read_csv("dataset2_all_chromosomes_exbor2.tsv", sep="\t")

### Filter the data
df1 = df1[df1['TE_dis'] <= 100000]
df2 = df2[df2['TE_dis'] <= 100000]
df3 = df3[df3['TE_dis'] <= 100000]
df4 = df4[df4['TE_dis'] <= 100000]

###top20_TEsuperfamilies = ["SINE/Alu", "LINE/L1", "SINE/MIR", "LINE/L2", "LTR/ERVL-MaLR", "DNA/hAT-Charlie", "LTR/ERV1", "LTR/ERVL", "DNA/TcMar-Tigger", "LINE/CR1", "DNA/hAT-Tip100", "DNA/hAT-Blackjack", "LTR/Gypsy", "DNA/TcMar-Mariner", "LINE/RTE-X", "LTR/ERVK", "LINE/RTE-BovB", "DNA/hAT", "DNA/TcMar-Tc2", "LTR/Gypsy?"]

### KDE plot
### Plotting KDE for all datasets with different colors
plt.figure(figsize=(10, 6))
sns.kdeplot(df1['TE_dis'], bw_adjust=0.1, label='Exbor1 Dataset1', color='blue')
sns.kdeplot(df2['TE_dis'], bw_adjust=0.1, label='Exbor2 Dataset1', color='red')
sns.kdeplot(df3['TE_dis'], bw_adjust=0.1, label='Exbor1 Dataset2', color='green')
sns.kdeplot(df4['TE_dis'], bw_adjust=0.1, label='Exbor2 Dataset2', color='purple')

plt.xlabel('Exbor')
plt.ylabel('Density')
plt.title('KDE: Merged datasets, all chr, both exbors: bw_adjust=0.1')
plt.grid(True)

### Set the x-axis to show ticks at intervals of 10000
plt.xticks(range(0, 100001, 10000))
plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))

### Add a legend
plt.legend()

### Show the plot
###plt.show()

# Save the plot as a PNG file
plt.savefig("KDE_plot.png", dpi=300)
plt.close()




### Broken-line plot
bins = np.arange(0, 100001, 400)

# Calculate counts for each bin
counts1, _ = np.histogram(df1['TE_dis'], bins=bins)
counts2, _ = np.histogram(df2['TE_dis'], bins=bins)
counts3, _ = np.histogram(df3['TE_dis'], bins=bins)
counts4, _ = np.histogram(df4['TE_dis'], bins=bins)

# Plotting broken line plot for all datasets with different colors
plt.figure(figsize=(10, 6))
plt.plot(bins[:-1], counts1, label='Exbor1 Dataset1', color='blue', marker='o', markersize=4, linestyle='--')
plt.plot(bins[:-1], counts2, label='Exbor2 Dataset1', color='red', marker='o', markersize=4, linestyle='--')
plt.plot(bins[:-1], counts3, label='Exbor1 Dataset2', color='green', marker='o', markersize=4, linestyle='--')
plt.plot(bins[:-1], counts4, label='Exbor2 Dataset2', color='purple', marker='o', markersize=4, linestyle='--')

plt.xlabel('Exbor')
plt.ylabel('Counts')
plt.title('Broken lines: Merged datasets, all chr, both exbors, 400nt bin size')
plt.grid(True)

# Set the x-axis to show ticks at intervals of 10000
plt.xticks(range(0, 100001, 10000))
plt.legend()
###plt.show()

# Save the plot as a PNG file
plt.savefig("Broken_lines_plot_400nt.png", dpi=300)
plt.close()