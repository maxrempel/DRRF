import pandas as pd
import re
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

def calculate_distance(row1, row2):
    position1 = row1['POS']
    cigar1 = row1['CIGAR']
    flag1 = row1['FLAG']
    position2 = row2['POS']
    cigar2 = row2['CIGAR']
    flag2 = row2['FLAG']
    dist = abs(position1 - position2)

    return flag1, position1, cigar1, flag2, position2, cigar2, dist

def print_rows_with_read_id(df, read_id, folder):
    rows_with_read_id = df.query("QNAME == @read_id")
    num_rows = len(rows_with_read_id)
    if num_rows < 2:
        # print(f"Only one row found for read ID {read_id}")
        return
    for i in range(num_rows - 1):
        with open(f"{folder}/splitted_long_dist.tsv", 'a') as f:
            FLAG1, pos1, cig1, FLAG2, pos2, cig2, distance = calculate_distance(rows_with_read_id.iloc[i], rows_with_read_id.iloc[i + 1])
            if distance >= 300:
                f.write(f"{read_id}\t{FLAG1}\t{pos1}\t{cig1}\t{FLAG2}\t{pos2}\t{cig2}\t{distance}\n")


def process_dist(folder_chr, chromosome):
    column_names = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR"]

    df = pd.read_csv(f"{folder_chr}/{chromosome}_bwa_band_no-mismatches.sam", sep='\t', usecols=range(6), names=column_names)
    df = df[df['CIGAR'].astype(str).str.contains(r'\d+M\d+S') | df['CIGAR'].astype(str).str.contains(r'\d+M\d+H')]

    # Read the file containing read IDs to search
    with open(f"{folder_chr}/splitted_IDs.lst", 'r') as f:
        read_ids_to_search = [line.strip() for line in f]

    if os.path.exists(f"{folder_chr}/splitted_long_dist.tsv"):
        os.remove(f"{folder_chr}/splitted_long_dist.tsv")
    subprocess.run('echo "read_id\tflag1\tpos1\tcigar1\tflag2\tpos2\tcigar2\tdistance" > {}/splitted_long_dist.tsv'.format(folder_chr), shell=True)

    # Loop through the read IDs and print the rows with each read ID
    for read_id in read_ids_to_search:
        with ThreadPoolExecutor() as executor:
            executor.submit(print_rows_with_read_id, df, read_id, folder_chr)










#