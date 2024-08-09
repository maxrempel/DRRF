import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.ticker import MaxNLocator
import re
from scipy.stats import gaussian_kde

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

top20_TEsuperfamilies = ["SINE/Alu", "LINE/L1", "SINE/MIR", "LINE/L2", "LTR/ERVL-MaLR", "DNA/hAT-Charlie", "LTR/ERV1", "LTR/ERVL", "DNA/TcMar-Tigger", "LINE/CR1"]

for TE_family in top20_TEsuperfamilies:
    out_TE_family = re.sub(r'/', '-', TE_family)

    # Assuming your dataframe is named df
    df1_spcTE = df1[df1['TE'] == TE_family]
    df2_spcTE = df2[df2['TE'] == TE_family]

    df3_spcTE = df3[df3['TE'] == TE_family]
    df4_spcTE = df4[df4['TE'] == TE_family]

    ### KDE plot
    ### Plotting KDE for all datasets with different colors
    # plt.figure(figsize=(10, 6))
    # sns.kdeplot(df1_spcTE['TE_dis'], bw_adjust=0.1, label='Exbor1 Dataset1', color='blue')
    # sns.kdeplot(df2_spcTE['TE_dis'], bw_adjust=0.1, label='Exbor2 Dataset1', color='red')
    # sns.kdeplot(df3_spcTE['TE_dis'], bw_adjust=0.1, label='Exbor1 Dataset2', color='green')
    # sns.kdeplot(df4_spcTE['TE_dis'], bw_adjust=0.1, label='Exbor2 Dataset2', color='purple')
    #
    # plt.xlabel('Exbor')
    # plt.ylabel('Density')
    # plt.title(f'{TE_family} KDE: Merged datasets, all chr, both exbors: bw_adjust=0.1')
    # plt.grid(True)

##################

    #modified version
    # Define a function to calculate and scale the KDE
    def compute_kde(data, bw_adjust, scale_factor=1):
        print(bw_adjust)
        kde = gaussian_kde(data, bw_method=0.025)
        x = np.linspace(min(data), max(data), 1000)
        y = kde(x) * scale_factor
        return x, y


    # Compute KDE for df1_spcTE and df2_spcTE with scaling
    x1, y1_scaled = compute_kde(df1_spcTE['TE_dis'], bw_adjust=0.1, scale_factor=1.5)
    x2, y2_scaled = compute_kde(df2_spcTE['TE_dis'], bw_adjust=0.1, scale_factor=1.5)

    # Compute KDE for df3_spcTE and df4_spcTE without scaling
    x3, y3 = compute_kde(df3_spcTE['TE_dis'], bw_adjust=0.1)
    x4, y4 = compute_kde(df4_spcTE['TE_dis'], bw_adjust=0.1)

    # Plot the scaled KDEs for df1_spcTE and df2_spcTE
    plt.figure(figsize=(10, 6))
    plt.plot(x1, y1_scaled, label='Exbor1 Dataset1 (scaled)', color='blue')
    plt.plot(x2, y2_scaled, label='Exbor2 Dataset1 (scaled)', color='red')

    # Plot the KDEs for df3_spcTE and df4_spcTE
    plt.plot(x3, y3, label='Exbor1 Dataset2', color='green')
    plt.plot(x4, y4, label='Exbor2 Dataset2', color='purple')

    ### Set the x-axis to show ticks at intervals of 10000
    plt.xticks(range(0, 100001, 10000))
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))

    ### Add a legend
    plt.legend()
    plt.xlabel('Exbor')
    plt.ylabel('Density')
    plt.title(f'{TE_family} KDE: Merged datasets, all chr, both exbors: bw_adjust=0.1')
    plt.grid(True)

    ### Show the plot
    ###plt.show()

    # Save the plot as a PNG file
    plt.savefig(f"KDE_{out_TE_family}_plot.png", dpi=300)
    plt.close()




    ### Broken-line plot
    breaking_range = 2000
    bins = np.arange(0, 100001, breaking_range)

    # Calculate counts for each bin
    counts1, _ = np.histogram(df1_spcTE['TE_dis'], bins=bins)
    counts2, _ = np.histogram(df2_spcTE['TE_dis'], bins=bins)
    counts3, _ = np.histogram(df3_spcTE['TE_dis'], bins=bins)
    counts4, _ = np.histogram(df4_spcTE['TE_dis'], bins=bins)

    # Plotting broken line plot for all datasets with different colors
    plt.figure(figsize=(10, 6))
    plt.plot(bins[:-1], counts1, label='Exbor1 Dataset1', color='blue', marker='o', markersize=4, linestyle='--')
    plt.plot(bins[:-1], counts2, label='Exbor2 Dataset1', color='red', marker='o', markersize=4, linestyle='--')
    plt.plot(bins[:-1], counts3, label='Exbor1 Dataset2', color='green', marker='o', markersize=4, linestyle='--')
    plt.plot(bins[:-1], counts4, label='Exbor2 Dataset2', color='purple', marker='o', markersize=4, linestyle='--')

    plt.xlabel('Exbor')
    plt.ylabel('Counts')
    plt.title(f'{TE_family} Broken lines: Merged datasets, all chr, both exbors, {breaking_range}nt bin size')
    plt.grid(True)

    # Set the x-axis to show ticks at intervals of 10000
    plt.xticks(range(0, 100001, 10000))
    plt.legend()
    ###plt.show()

    # Save the plot as a PNG file
    plt.savefig(f'Broken_lines_{out_TE_family}_plot_{breaking_range}nt.png', dpi=300)
    plt.close()