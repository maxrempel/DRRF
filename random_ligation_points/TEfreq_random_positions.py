import pandas as pd
pd.set_option('display.max_columns', None)
from tqdm import tqdm
import re
from io import StringIO
import pybedtools
from concurrent.futures import ThreadPoolExecutor, as_completed
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import os
import numpy as np

def process_row(row, chr, TEbed):
    harbor1 = row[14]
    harbor2 = row[19]

    # cigar1 = int(arm_positions.iloc[0, 8])
    # cigar2 = int(arm_positions.iloc[0, 9])

    ##harbor1 = row[12]
    ##harbor2 = row[18]
    upstream_H1 = harbor1 - 100000
    if upstream_H1 < 0:
        upstream_H1 = 1
    downstream_H1 = harbor1 + 100000

    upstream_H2 = harbor2 - 100000
    if upstream_H2 < 0:
        upstream_H2 = 1
    downstream_H2 = harbor2 + 100000

    harbors_bed = (f"{chr}\t{upstream_H1}\t{downstream_H1}\tharbor1\n"
                   f"{chr}\t{upstream_H2}\t{downstream_H2}\tharbor2")
    harbors_pb = pybedtools.BedTool(harbors_bed, from_string=True)
    intersection = harbors_pb.intersect(TEbed, wa=True, wb=True, nonamecheck=True)

    if intersection:
        intersection_str = str(intersection)
        intersection_df = pd.read_csv(StringIO(intersection_str), sep='\t', header=None)
        intersection_df[5] = pd.to_numeric(intersection_df[5])
        intersection_df[6] = pd.to_numeric(intersection_df[6])
        intersection_df[10] = ((intersection_df[6] - intersection_df[5]) / 2).round().astype(int)
        intersection_df[11] = intersection_df[10] + intersection_df[5]
        if harbor1 < harbor2:
            intersection_df[12] = abs(intersection_df[11] - intersection_df[1])
        elif harbor1 > harbor2:
            intersection_df[12] = abs(intersection_df[11] - intersection_df[2])
        ### intersection_df.loc[intersection_df[3] == 'harbor1', 12] = ((intersection_df.loc[intersection_df[3] == 'harbor1', 11] - intersection_df.loc[intersection_df[3] == 'harbor1', 1]) - cigar1)
        ### intersection_df.loc[intersection_df[3] == 'harbor2', 12] = (intersection_df.loc[intersection_df[3] == 'harbor2', 11] - intersection_df.loc[intersection_df[3] == 'harbor2', 1])
        ###print(intersection_df)

        return intersection_df

    return None

def get_chr_size(chr_n):
    # Read the TSV file into a DataFrame
    df = pd.read_csv("chr_size.tsv", sep='\t', header=None, names=['Column1', 'Column2'])

    # Filter the DataFrame for the row where Column1 matches the var
    value_row = df[df['Column1'] == chr_n]

    # If the row is found, return the value from Column2
    if not value_row.empty:
        return value_row['Column2'].values[0]
    else:
        return None

def generate_random_values(chr_size):
    end_LG1 = np.random.randint(0, chr_size - 300)
    end_harbor2 = np.random.randint(end_LG1 + 300, chr_size)
    return end_LG1, end_harbor2

print("Performing analysis chromosome per chromosome...")

dataset_number = "1"
chr_list = ['chr' + str(i) for i in range(1, 23)] + ['chrX']

for chr in chr_list:
    if not os.path.exists(str(f"{chr}_out")):
        os.makedirs(f"{chr}_out")
    LP_output = pd.read_csv(f"dataset{dataset_number}_outputs/{chr}_MASKED_LigP_finder_output_opposite.tsv", sep='\t')
    chim_reads = pd.read_csv(f"dataset{dataset_number}_outputs/{chr}_chimeric_reads_rmDUP.tsv", sep='\t', header=None)

    #### Random LPS
    chr_size = get_chr_size(chr)

    for i in range(len(LP_output)):
        LP_output.at[i, 'end_LG1'], LP_output.at[i, 'end_harbor2'] = generate_random_values(chr_size)

    ####LP_output.to_csv("test_random.csv", index = False)
    ##Cigar 1
    chim_reads[8] = chim_reads[4]
    chim_reads[8] = chim_reads[8].apply(
        lambda x: int(re.search(r'(\d+)M', x).group(1)) if pd.notnull(x) and re.search(r'(\d+)M', x) else None)

    ##Cigar 2
    chim_reads[9] = chim_reads[7]
    chim_reads[9] = chim_reads[9].apply(
        lambda x: int(re.search(r'(\d+)M', x).group(1)) if pd.notnull(x) and re.search(r'(\d+)M', x) else None)

    ##Harbor1 position
    chim_reads[10] = pd.to_numeric(chim_reads[8] + chim_reads[3])
    ##Harbor2 position
    chim_reads[11] = pd.to_numeric(chim_reads[9] + chim_reads[6] - 1)

    chim_reads.to_csv(f"{chr}_out/{chr}_chim_reads_harbor_match.tsv", header=True, index=False, sep='\t')

    TE_fam_pybed = pybedtools.BedTool(f"hg38_only_TEs.bed")

    ### Create a ThreadPoolExecutor
    num_workers = 6  # Adjust the number of workers based on your system's capability
    results = []
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(process_row, row, chr, TE_fam_pybed) for row in LP_output.itertuples(index=False)]

        for future in tqdm(as_completed(futures), total=len(futures)):
            result = future.result()
            if result is not None:
                results.append(result)

    ### Concatenate all DataFrames and write once
    if results:
        intersection_df = pd.concat(results, ignore_index=True)

    intersection_df.columns = ['H_chr', 'H_start', 'H_end', 'Harbor', 'TE_chr', 'TE_start', 'TE_end',
                               'TE', 'score', 'TE_strand', 'TE_half', 'TE_middle_pos', 'TE_dis']
    intersection_df = intersection_df.drop_duplicates()

    harbor1 = intersection_df[intersection_df['Harbor'] == 'harbor1']

    ### Calculate the 30th percentile of the score column
    score_threshold = harbor1['score'].quantile(0.3)
    ### Filter the DataFrame to exclude rows with scores below the threshold
    harbor1 = harbor1[harbor1['score'] > score_threshold]

    harbor2 = intersection_df[intersection_df['Harbor'] == 'harbor2']
    ### Calculate the 30th percentile of the score column
    score_threshold = harbor2['score'].quantile(0.3)
    ### Filter the DataFrame to exclude rows with scores below the threshold
    harbor2 = harbor2[harbor2['score'] > score_threshold]

    if os.path.exists(str(f"{chr}_out/{chr}_harbor1.tsv")):
        os.remove(str(f"{chr}_out/{chr}_harbor1.tsv"))
    if os.path.exists(str(f"{chr}_out/{chr}_harbor2.tsv")):
        os.remove(str(f"{chr}_out/{chr}_harbor2.tsv"))

    harbor1.to_csv(f"{chr}_out/{chr}_harbor1.tsv", header=True, index=False, sep='\t')
    harbor2.to_csv(f"{chr}_out/{chr}_harbor2.tsv", header=True, index=False, sep='\t')


### Merge all chromosomes per TE superfamily - per exbor
all_chr_harbor1 = []
all_chr_harbor2 = []
for chr_id in chr_list:
    file_path = f"{chr_id}_out/{chr_id}_harbor1.tsv"
    all_chr_harbor1.append(pd.read_csv(file_path, sep="\t"))
    merged_chromosomes_harbor1 = pd.concat(all_chr_harbor1, ignore_index=True)

    file_path = f"{chr_id}_out/{chr_id}_harbor2.tsv"
    all_chr_harbor2.append(pd.read_csv(file_path, sep="\t"))
    merged_chromosomes_harbor2 = pd.concat(all_chr_harbor2, ignore_index=True)

merged_chromosomes_harbor1.to_csv(f"dataset{dataset_number}_all_chromosomes_exbor1.tsv", sep="\t", index=False, header=True)
merged_chromosomes_harbor2.to_csv(f"dataset{dataset_number}_all_chromosomes_exbor2.tsv", sep="\t", index=False, header=True)