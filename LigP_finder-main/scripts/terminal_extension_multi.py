import warnings
from Bio import BiopythonWarning
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=BiopythonWarning)
    from Bio import pairwise2
from Bio.Seq import Seq
import tempfile
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
import pandas as pd
from __main__ import out_folder
import os
import argparse
from __main__ import folder_chr
from __main__ import check_file
from __main__ import genome_file

#genome_file = f"{genome}"
pattern = r'(\d+)M'

def get_sequence_from_fastq(read_id, fq_file):
    from __main__ import out_folder
    lines = fq_file.readlines()
    for i in range(0, len(lines), 4):
        header = lines[i].strip()
        sequence = lines[i+1].strip()
        if header.split()[0][1:] == read_id:
            return sequence
    return None

def reverse_complement(seq):
    return seq.reverse_complement()

def count_lines(string):
    # Split the string by newline characters and count the number of resulting substrings
    return len(string.split('\n'))

def process_record(row, sequence, fastq):
    ID = row[0]
    # if ID == "SRR12625672.15760902":
    ##### Terminal fragment extension
    term_START = row[6]
    term_long = row[7]
    term_long = re.search(pattern, term_long)
    term_long = int(term_long.group(1))
    term_END = (term_START + int(term_long)) - 1
    upstream_term = term_START - 50
    downstream_term = term_END + 50
    term_frag = sequence[term_START - 1:term_END]
    ext_term = sequence[upstream_term - 1:downstream_term]

    ##### Internal fragment extension
    int_START = row[3]
    int_long = row[4]
    int_long = re.search(pattern, int_long)
    int_long = int(int_long.group(1))
    int_END = (int_START + int(int_long)) - 1
    upstream_int = int_START - 50
    downstream_int = int_END + 50
    int_frag = sequence[int_START - 1:int_END]
    ext_int = sequence[upstream_int - 1:downstream_int]

    int_frag_len = len(int_frag)
    term_frag_len = len(term_frag)


    if int_frag_len + term_frag_len >= 149 and int_frag_len + term_frag_len <= 151:
        qcov = abs((term_frag_len * 100) / 150)
        read_sequence = get_sequence_from_fastq(ID, fastq)

        with tempfile.NamedTemporaryFile(prefix='ext_frag_', suffix='.fa', delete=True) as ext_frag_file, \
                tempfile.NamedTemporaryFile(prefix='int_frag_', suffix='.fa', delete=True) as int_frag_file, \
                tempfile.NamedTemporaryFile(prefix='chim_read_', suffix='.fa', delete=True) as chimeric_read_file:

            #### Write sequence data to the temporary files
            with open(ext_frag_file.name, 'w') as f:
                f.write(f">extended_frag\n{ext_term}\n")

            with open(int_frag_file.name, 'w') as f2:
                f2.write(f">extended_frag\n{ext_int}\n")

            with open(chimeric_read_file.name, 'w') as f3:
                f3.write(f">chim_read\n{read_sequence}\n")

            command = 'blastn -task blastn-short -query {} -subject {} -perc_identity 100 -outfmt "6 qseqid sseqid pident sstart send evalue sseq" -ungapped -qcov_hsp_perc {}'.format(
                chimeric_read_file.name, ext_frag_file.name, qcov)
            result_terminal = subprocess.run(command, shell=True, capture_output=True, text=True)

            command = 'blastn -task blastn-short -query {} -subject {} -perc_identity 100 -outfmt "6 qseqid sseqid pident sstart send evalue sseq" -ungapped -qcov_hsp_perc {}'.format(
                chimeric_read_file.name, int_frag_file.name, qcov)
            result_internal = subprocess.run(command, shell=True, capture_output=True, text=True)

            if result_terminal.returncode == 0 and result_internal.returncode == 0:
                terminal_out = result_terminal.stdout
                lines = terminal_out.strip().split('\n')
                last_column_values = [line.split('\t')[-1] for line in lines]
                terminal_seq = '\t'.join(last_column_values)
                terminal_len = len(terminal_seq)

                if terminal_len > 0:
                    n_matches_terminal = count_lines(terminal_seq)

                    internal_out = result_internal.stdout
                    lines = internal_out.strip().split('\n')
                    last_column_values = [line.split('\t')[-1] for line in lines]
                    internal_seq = '\t'.join(last_column_values)
                    internal_len = len(internal_seq)
                    n_matches_internal = count_lines(internal_seq)


                    if n_matches_terminal == 1 and n_matches_internal == 1:
                        if terminal_len == term_frag_len and internal_len == int_frag_len:

                            with open(f"{folder_chr}/result_ligation_points.lst", 'a') as out:
                                out.write(f"{ID}\n")


if __name__ == 'terminal_extension_multi':
    if check_file(f"{folder_chr}/result_ligation_points.lst") == True:
        os.remove(f"{folder_chr}/result_ligation_points.lst")
    chimeric_reads = pd.read_csv(f"{folder_chr}/chimeric_reads.tsv", sep='\t')
    chimeric_reads = chimeric_reads.drop_duplicates(subset=chimeric_reads.columns[1:7], keep='first')

    with open(f"{folder_chr}/SGreads_merged.fastq", 'r') as fastq_file:
        with open(genome_file, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        with ThreadPoolExecutor() as executor:
            for row in chimeric_reads.itertuples(index=False):
                for record in records:
                    sequence = record.seq
                    executor.submit(process_record, row, sequence, fastq_file)