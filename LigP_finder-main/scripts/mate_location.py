import pandas as pd
import subprocess
import os


def process(folder_chr, chromosome):
    print(f"Folder_chr = {folder_chr}")
    print(f"chromosome = {chromosome}")
    SGmates = pd.read_csv(f"{folder_chr}/{chromosome}_SGmate.sam", sep='\t', header=None, usecols=[0, 3])
    SGmates.columns = ['ID_SRR', 'pos_SRR']
    print(SGmates)
    splitted_reads = pd.read_csv(f"{folder_chr}/splitted_long_dist.tsv", sep='\t')
    print(splitted_reads)
    # Create an empty list to store rows from df2
    rows_with_same_value = []

    if os.path.exists(f"{folder_chr}/chimeric_reads.tsv"):
        print("chimeric_read.tsv exists! remove it")
        os.remove(f"{folder_chr}/chimeric_reads.tsv")

    subprocess.run('echo "read_id\tcontinuous_read\tflag1\tinternal_arm\tinternal_cigar\tflag2\tterminal_arm\tterminal_cigar" > {}/chimeric_reads.tsv'.format(folder_chr), shell=True)

    print("creating chimeric_reads.tsv")
    for row in splitted_reads.itertuples(index=False):
        # Extract the value from the first column of the current row
        read_ID = row[0]
        # Filter rows from df2 based on the value of the first column
        position_mate = SGmates.query("ID_SRR == @read_ID")
        position_mate = position_mate.iloc[:, 1].values[0]
        pos_mate_plus = abs(int(position_mate) + 500)
        pos_mate_minus = abs(int(position_mate) - 500)

        flag1 = row[1]
        pos_frag1 = row[2]
        cigar_frag1 = row[3]

        flag2 = row[4]
        pos_frag2 = row[5]
        cigar_frag2 = row[6]

        with open(f"{folder_chr}/chimeric_reads.tsv", 'a') as f:
            if pos_frag1 < pos_mate_plus and pos_frag1 > pos_mate_minus:
                f.write(f"{read_ID}\t{position_mate}\t{flag1}\t{pos_frag1}\t{cigar_frag1}\t{flag2}\t{pos_frag2}\t{cigar_frag2}\n")
            elif pos_frag2 < pos_mate_plus and pos_frag2 > pos_mate_minus:
                f.write(f"{read_ID}\t{position_mate}\t{flag1}\t{pos_frag1}\t{cigar_frag1}\t{flag2}\t{pos_frag2}\t{cigar_frag2}\n")




#