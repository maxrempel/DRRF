import subprocess
import pandas as pd
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

pattern = r'(\d+)M'
genome_file = "hg38_chr20_unmasked.fa"

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def calc_ligdist(result_blast, start_posLG1, end_posLG1, start_posLG2, end_posLG2, out_df):
    # Convert tuple to pandas Series
    series = pd.Series(result_blast)
    # Create a DataFrame from the Series
    result_blast = pd.DataFrame(series).T  # Transpose to make it row-wise
    # Rename columns
    result_blast.columns = ['query_ID', 'subject_ID', 'length', 'perc_id', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'q_seq', 's_seq']


    #Fixing reverse matches: start position always smaller than end position
    q_start = result_blast['q_start'].item()
    q_end = result_blast['q_end'].item()
    s_start = result_blast['s_start'].item()
    s_end = result_blast['s_end'].item()
    if q_start > q_end:
        result_blast['q_start'] = q_end
        result_blast['q_end'] = q_start
    if s_start > s_end:
        result_blast['s_start'] = s_end
        result_blast['s_end'] = s_start

    ########## LP1
    result_blast['start_harbor'] = start_posLG1
    result_blast['end_harbor'] = end_posLG1

    result_blast['end_LG1'] = result_blast['end_harbor'] + 1
    result_blast['start_homolog1'] = result_blast['start_harbor'] + result_blast['q_start'].item()
    result_blast['end_homolog1'] = result_blast['start_harbor'] + result_blast['q_end'].item()

    ##distance LP1
    result_blast['LG1_distance'] = result_blast['end_harbor'] - result_blast['end_homolog1']


    ########## LP2
    result_blast['start_harbor2'] = start_posLG2
    result_blast['end_harbor2'] = end_posLG2
    result_blast['end_LG2'] = result_blast['start_harbor2'] + 1
    result_blast['start_homolog2'] = result_blast['start_harbor2'] + result_blast['s_start'].item()
    result_blast['end_homolog2'] = result_blast['start_harbor2'] + result_blast['s_end'].item()


    ##distance LG2
    #IMPORTANT IF REVERSE
    result_blast['LG2_distance'] = result_blast['end_harbor2'] - result_blast['end_homolog2']

    ### IF PLUS
    # result_blast['LG2_distance'] = result_blast['start_homolog2'] - result_blast['start_harbor2']

    ###Calculate differential difference
    Ligdist = abs(result_blast['LG2_distance'] - result_blast['LG1_distance']).item()
    if Ligdist < 6:
        # upstream_HARBOR = reverse_complement(upstream_HARBOR)
        result_blast['Ligdis'] = Ligdist
        if out_df is None:
            out_df = result_blast
        else:
            out_df = out_df.append(row, ignore_index=True)

    if out_df is not None:
        return out_df

def homologs_blast(row, sequence):
    ### Internal arm
    cigar_int = row[2]
    cigar_int = re.search(pattern, cigar_int)
    cigar_int = int(cigar_int.group(1))
    start_pos_int = row[1]
    # if start_pos_int == 36332118:

    ### Terminal arm
    cigar_term = row[5]
    cigar_term = re.search(pattern, cigar_term)
    cigar_term = int(cigar_term.group(1))
    start_pos_term = row[4]

    if start_pos_int > start_pos_term:
        start_pos_int, start_pos_term = start_pos_term, start_pos_int

    # if start_pos_int == 49641255:
    end_pos_int = (start_pos_int + cigar_int) - 1
    start_upstream = end_pos_int - 120
    upstream_LG_primary = sequence[start_upstream - 1:end_pos_int]
    end_pos_term = (start_pos_term + cigar_term) - 1
    end_downstream = end_pos_term + 120
    downstream_LG_primary = sequence[end_pos_term - 1:end_downstream]


    #### Purine code here
    upstream_LG = re.sub(r'a|A', 'G', str(upstream_LG_primary))
    upstream_LG = re.sub(r't|T', 'C', str(upstream_LG))

    downstream_LG = re.sub(r'a|A', 'G', str(downstream_LG_primary))
    downstream_LG = re.sub(r't|T', 'C', str(downstream_LG))

    ######### IMPORTANT STEP, CHANGES EVERYTHING
    ### Reverse complement downstream harbor - IF REVERSE
    downstream_LG = reverse_complement(downstream_LG)

    with tempfile.NamedTemporaryFile(prefix='ext_frag_', suffix='.fa', delete=True) as upstream_LG_f, \
            tempfile.NamedTemporaryFile(prefix='int_frag_', suffix='.fa', delete=True) as downstream_LG_f:

        with open(upstream_LG_f.name, 'w') as f:
            f.write(f">upstream_LG\n{upstream_LG}\n")
        with open(downstream_LG_f.name, 'w') as f2:
            f2.write(f">downstream_LG\n{downstream_LG}\n")

        command = 'blastn -task blastn-short -query {} -subject {} -perc_identity 80 -qcov_hsp_perc 27 -ungapped -outfmt "6 qseqid sseqid length pident qstart qend sstart send evalue qseq sseq"'.format(
        upstream_LG_f.name, downstream_LG_f.name)
        res_blast = subprocess.run(command, shell=True, capture_output=True, text=True)
        if res_blast.returncode == 0:
            # Convert the output string into a file-like object
            output_string = res_blast.stdout
            output_file_like_object = StringIO(output_string)
            # Read the file-like object using pandas
            result_blast_df = pd.read_csv(output_file_like_object, sep='\t', header=None)
            result_blast_df.columns = ['query_ID', 'subject_ID', 'length', 'perc_id', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'q_seq', 's_seq']
            appended_df = None
            for row in result_blast_df.itertuples(index=False):
                satisfied_homologs = calc_ligdist(row, start_upstream, end_pos_int, end_pos_term, end_downstream, appended_df)

            if satisfied_homologs is not None:
                first_row_values = satisfied_homologs.loc[0, ['end_harbor', 'start_harbor2']]
                # Concatenating the values with an underscore
                out_file = '_'.join(first_row_values.astype(str))

                with tempfile.NamedTemporaryFile(prefix='blast_out', suffix='.tsv', delete=True) as blast_table:
                    satisfied_homologs.to_csv(f"{blast_table.name}", header=True, index=False, sep='\t')
                    subprocess.run(["Rscript", "--vanilla", "dotplot.R", blast_table.name,
                                    out_file])
                # Repeat the string value for each row in the DataFrame
                new_column = pd.DataFrame({'LigP_ID': [f"LigP_{out_file}"] * len(satisfied_homologs)})

                # Concatenate the new DataFrame with the original DataFrame, placing the new column at the beginning
                satisfied_homologs = pd.concat([new_column, satisfied_homologs], axis=1)


                print('Dataframe:\n', satisfied_homologs)
                satisfied_homologs.to_csv('LigP_finder_output.tsv', mode='a', header=False, index=False, sep = '\t')





            # for row in df.itertuples(index=False):
            #     plot_LG1(row, start_upstream, end_pos_int, end_pos_term, end_downstream, start_pos_int, upstream_LG, downstream_LG, upstream_LG_primary, downstream_LG_primary)


if __name__ == '__main__':
    chimeric_reads = pd.read_csv('result_ligation_points_rmDUP.tsv', sep='\t', header=None)
    with open(genome_file, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    with ThreadPoolExecutor() as executor:
        for row in chimeric_reads.itertuples(index=False):
            for record in records:
                sequence = record.seq
                executor.submit(homologs_blast, row, sequence)
    if os.path.exists('LigP_finder_output.tsv'):
        df = pd.read_csv('LigP_finder_output.tsv', sep='\t', header=None)
        df.columns = ['LigP_ID', 'query_ID', 'subject_ID', 'length', 'perc_id', 'q_start', 'q_end' , 's_start',
                      's_end', 'evalue', 'q_seq', 's_seq', 'start_harbor', 'end_harbor', 'end_LG1',
                      'start_homolog1', 'end_homolog1',	'LG1_distance',	'start_harbor2', 'end_harbor2',
                      'end_LG2', 'start_homolog2', 'end_homolog2', 'LG2_distance', 'Ligdis']
        df.to_csv('LigP_finder_output.tsv', header=True, index=False, sep='\t')





#