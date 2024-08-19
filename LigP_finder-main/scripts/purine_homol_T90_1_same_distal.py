import subprocess
import pandas as pd
from tqdm import tqdm
import tempfile
import os
import re
import random
import string
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from io import StringIO
import warnings
from Bio import BiopythonWarning
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=BiopythonWarning)

from __main__ import spotcheck_continous

if __name__ == 'purine_homol_T90_1_same_distal':
    def process_same(folder_chr, chromosome):
        chimeric_reads = pd.read_csv(f"{folder_chr}/chimeric_reads.tsv", sep="\t")
        chimeric_reads = chimeric_reads.drop_duplicates(subset=chimeric_reads.columns[1:7], keep='first')

        chimeric_reads.to_csv(f"{folder_chr}/chimeric_reads_rmDUP.tsv", header=False, index=False, sep="\t")
        chimeric_reads = pd.read_csv(f"{folder_chr}/chimeric_reads_rmDUP.tsv", sep='\t')

        chimeric_reads.columns = ['read_id', 'continuous_read', 'flag1', 'internal_arm', 'internal_cigar', 'flag2',
                                  'terminal_arm', 'terminal_cigar']
        chimeric_reads = chimeric_reads[chimeric_reads['internal_cigar'].str.match(r'\d+M\d+[SH]')]
        chimeric_reads = chimeric_reads[chimeric_reads['terminal_cigar'].str.match(r'\d+M\d+[SH]')]

        chimeric_reads = chimeric_reads[
            chimeric_reads['internal_cigar'].astype(str).str.contains(r'\d+M\d+S') & chimeric_reads[
                'terminal_cigar'].astype(str).str.contains(r'\d+M\d+H')]

        #### Invert values where terminal_arm < continuous_read
        mask = chimeric_reads['terminal_arm'] < chimeric_reads['internal_arm']
        chimeric_reads.loc[mask, ['internal_arm', 'terminal_arm', 'internal_cigar', 'terminal_cigar']] = \
        chimeric_reads.loc[mask, ['terminal_arm', 'internal_arm', 'terminal_cigar', 'internal_cigar']].values
        # chimeric_reads.loc[mask, ['internal_arm', 'terminal_arm']] = chimeric_reads.loc[
        #     mask, ['terminal_arm', 'internal_arm']].values


        chimeric_reads['distance'] = abs(chimeric_reads['terminal_arm'] - chimeric_reads['continuous_read'])
        chimeric_reads = chimeric_reads[chimeric_reads['distance'] > 450]
        ##chimeric_reads = chimeric_reads[chimeric_reads['continuous_read'] < 145138636]

        chimeric_reads['continuous_read'] = abs(chimeric_reads['continuous_read'] + 150)
        chimeric_reads['upstream_cont'] = abs(chimeric_reads['continuous_read'] - 299)
        chimeric_reads['downstream_LP'] = abs(chimeric_reads['terminal_arm'] - 299)
        chimeric_reads.to_csv(f"{folder_chr}/chimeric_reads_distance_distal.tsv", header=True, index=False, sep="\t")


        with tqdm(total=len(chimeric_reads)) as pbar:
            with ThreadPoolExecutor() as executor:
                for row in chimeric_reads.itertuples(index=False):
                    executor.submit(spotcheck_continous, row, str("opposite"), pbar, folder_chr, chromosome)
            if os.path.exists(f"{folder_chr}/dotplots_opposite/distal_LigP_finder_output_opposite.tsv"):
                df = pd.read_csv(f"{folder_chr}/dotplots_opposite/distal_LigP_finder_output_opposite.tsv", sep='\t', header=None)
                df.columns = ['LigP_ID', 'query_ID', 'subject_ID', 'length', 'perc_id', 'q_start', 'q_end' , 's_start',
                              's_end', 'evalue', 'q_seq', 's_seq', 'start_harbor', 'end_harbor', 'end_LG1',
                              'start_homolog1', 'end_homolog1',	'LG1_distance',	'start_harbor2', 'end_harbor2',
                              'end_LG2', 'start_homolog2', 'end_homolog2', 'LG2_distance', 'Ligshift']
                df = df.drop_duplicates()
                df.to_csv(f"{folder_chr}/dotplots_opposite/distal_LigP_finder_output_opposite.tsv", header=True, index=False, sep='\t')





#
