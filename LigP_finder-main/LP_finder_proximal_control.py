import subprocess
import os
import sys
import re
import tempfile
import argparse
from io import StringIO
from Bio import SeqIO
import pandas as pd
import pysam

pd.set_option('display.max_columns', None)

sys.path.insert(1, 'scripts/')
out_folder = str('data/')

parser = argparse.ArgumentParser(description='This script has been designed to identify ligation points, extract harbors, and find homologs between them.')

required = parser.add_argument_group('Required arguments')
required.add_argument('--genome', help='Genome in fasta', required=True, type=str, metavar="")
required.add_argument('--read1', help='Fastq R1 from paired-end', required=True, type=str, metavar="")
required.add_argument('--read2', help='Fastq R2 from paired-end', required=True, type=str, metavar="")

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('--harbor_length', help='Harbor length in nucleotides. Default = 300', required=False, type=int, default=300, metavar="")
optional.add_argument('--cov', help='Minimum percentage length of the harbor to the homologs. Default = 10', required=False, type=float, default=10, metavar="")
optional.add_argument('--perc', help='Minimum homology percentage between homologs. Default = 80', required=False, type=int, default=80, metavar="")
optional.add_argument('--out', help='Folder name to make the output files. Default = output_prox_spotcheck', required=False, type=str, default="LP_finder_output", metavar="")
optional.add_argument('--dotplot', help='Activate the creation of 50 dotplots per orientation', required=False, action = 'store_true')

args = parser.parse_args()

genome_file = args.genome
out_folder = str(f"{args.out}/")
chr_number = re.sub(r'output_|/', '', str(out_folder))

def check_file(f):
    if os.path.exists(f):
        if os.stat(str(f)).st_size > 0:
            return True
        else:
            return False
    else:
        return False
if check_file(f"{args.out}") == False:
    os.mkdir(f"{args.out}")

def complementar(seq, op_harbor):
    if op_harbor:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '>>':'<<'}
    else:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    #return ''.join(complement.get(base, base) for base in reversed(seq))
    return ''.join(complement.get(base, base) for base in seq)

def complementar_rev(seq):
    complement = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def complementar_and_rev(seq):
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


    ### distance LP2
    result_blast['LG2_distance'] = result_blast['end_harbor2'] - result_blast['end_homolog2']

    ###Calculate differential difference
    Ligdist = abs(result_blast['LG2_distance'] - result_blast['LG1_distance']).item()
    # if Ligdist < 6:
    if Ligdist:
        # upstream_HARBOR = reverse_complement(upstream_HARBOR)
        result_blast['Ligdis'] = Ligdist
        if out_df is None:
            out_df = result_blast
        else:
            out_df = out_df.append(row, ignore_index=True)

    if out_df is not None:
        return out_df

def insert_brackets(harbor_seq, start, end, strandness, op_harbor):
    if not op_harbor:
        harbor_marked = harbor_seq[:start + 1] + ">>" + harbor_seq[start - 1:]
        harbor_marked = harbor_marked[:end + 4] + ">>" + harbor_marked[end + 4:]
    elif op_harbor:
        harbor_marked = harbor_seq[:start - 1] + "<<" + harbor_seq[start -1:]
        harbor_marked = harbor_marked[:end + 2] + "<<" + harbor_marked[end + 3:]

    return harbor_marked


def extract_sequence_by_id(fastq_file, id):
    # Open the FASTQ file and iterate over records
    with open(fastq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # Check if the record ID matches the specified ID
            if record.id == id:
                return str(record.seq)  # Return the sequence as a string

def primary2purine(seq):
    purine = re.sub(r'a|A', 'G', str(seq))
    purine = re.sub(r't|T', 'C', str(purine))
    return purine

def extract_sequence(start, end):
    fasta_file = f"data/hg38_chr8_unmasked.fa"
    faidx_file = fasta_file + ".fai"

    # Create a FastaFile object from indexed FASTA file
    fasta = pysam.FastaFile(fasta_file)

    # Check if chromosome exists in the index
    #if chromosome not in fasta.references:
    #    print(f"Chromosome {chromosome} not found in the file.")
    #    return ""
    # Extract sequence
    #if strand == '+':
    sequence = fasta.fetch(str("chr8"), start - 1, end).upper()
    #elif strand == '-':
    #    sequence = fasta.fetch(chromosome, start - 1, end).upper()
    #    sequence = complementar_and_rev(sequence)
    # Close the FastaFile object
    fasta.close()
    return sequence

def spotcheck_continous(row, strand, load_pbar, out_path, chr_n):
    load_pbar.update(1)
    read_ID = row[0]
    start_pos_continous_RAW = row[1]
    start_pos_term_RAW = row[6]
    ##if read_ID == "SRR12625672.10171203":
    #### Getting read sequence
    chim_reads = pd.read_csv(f"{out_path}/chimeric_reads_distance_SPOTCHECK.tsv", sep='\t')
    #### Get read split into the two positions
    read_ID = chim_reads[(chim_reads['continuous_read'] == start_pos_continous_RAW) & (chim_reads['terminal_arm'] == start_pos_term_RAW)].iloc[0]
    read_ID = read_ID['read_id']

    chim_read_sequence = extract_sequence_by_id(f"{out_path}/{chr_n}_SGreads_R1.fastq", read_ID)
    if chim_read_sequence is not None:
        cont_read_sequence = extract_sequence_by_id(f"{out_path}/all_chim_reads_R2.fastq", read_ID)
    elif chim_read_sequence is None:
        chim_read_sequence = extract_sequence_by_id(f"{out_path}/{chr_n}_SGreads_R2.fastq", read_ID)
        cont_read_sequence = extract_sequence_by_id(f"{out_path}/all_chim_reads_R1.fastq", read_ID)

    ####Get raw read sequence
    chim_read_sequence_rev = complementar(chim_read_sequence, False)
    cont_read_sequence_rev = complementar(cont_read_sequence, False)

    ### Terminal arm length
    cigar_term = row[7]
    cigar_term = re.search(r'(\d+)M', cigar_term)
    cigar_term = int(cigar_term.group(1))

    ### Internal arm length
    cigar_int = cigar_term
    ##print("cigar_int ", cigar_int)

    # if start_pos_int_RAW > start_pos_term_RAW:
    #     start_pos_int_RAW, start_pos_term_RAW = start_pos_term_RAW, start_pos_int_RAW
    #     cigar_int, cigar_term = cigar_term, cigar_int

    ### Get read pieces continous read
    cont_read_internal = cont_read_sequence[:cigar_int]
    cont_read_internal_rev = complementar(cont_read_internal, False)
    ##print(f"Piece 1:\n{cont_read_internal}")

    cont_read_terminal = cont_read_sequence[cigar_int:]
    cont_read_terminal_rev = complementar(cont_read_terminal, False)

    ### Get read pieces chimeric read
    chim_read_internal = chim_read_sequence[:cigar_term]
    chim_read_internal_rev = complementar(chim_read_internal, False)
    ##print(f"Piece 1:\n{chim_read_internal}")

    chim_read_terminal = chim_read_sequence[cigar_term:]
    chim_read_terminal_rev = complementar(chim_read_terminal, False)

    ### Test if continous read has 150nt long match
    with tempfile.NamedTemporaryFile(prefix='cont_read_', suffix='.fa', delete=True) as cont_read:
        with open(cont_read.name, 'w') as f:
            f.write(f">cont_read\n{cont_read_sequence}\n")
        command = 'blastn -task blastn-short -ungapped -query {} -subject {} -qcov_hsp_perc 100 -outfmt "6 qseqid sseqid length pident qstart qend sstart send evalue qseq sseq"'.format(
            cont_read.name, str(f"{out_path}/{chr_n}.fa"))
        cont_read_blast = subprocess.run(command, shell=True, capture_output=True, text=True)
        output_string = cont_read_blast.stdout
        if output_string.strip():
            output_file_like_object = StringIO(output_string)
            cont_read_blast_df = pd.read_csv(output_file_like_object, sep='\t', header=None)
            cont_read_blast_df.columns = ['query_ID', 'subject_ID', 'length', 'perc_id', 'q_start', 'q_end',
                                          's_start',
                                          's_end', 'evalue', 'q_seq', 's_seq']

            cont_length = cont_read_blast_df.loc[0, 'length']
            if cont_length == 150:
                ### Position piece 2 from continous read
                with tempfile.NamedTemporaryFile(prefix='piece2_', suffix='.fa', delete=True) as piece2_pos:
                    with open(piece2_pos.name, 'w') as f:
                        f.write(f">piece2\n{cont_read_terminal}\n")
                    command = 'blastn -task blastn-short -ungapped -query {} -subject {} -qcov_hsp_perc 100 -outfmt "6 qseqid sseqid length pident qstart qend sstart send evalue qseq sseq"'.format(
                        piece2_pos.name, str(f"{out_path}/{chr_n}.fa"))
                    piece2_blast = subprocess.run(command, shell=True, capture_output=True, text=True)
                    #### Convert the output string into a file-like object
                    output_string = piece2_blast.stdout
                    if output_string.strip():
                        output_file_like_object = StringIO(output_string)
                        piece2_blast_df = pd.read_csv(output_file_like_object, sep='\t', header=None)
                        piece2_blast_df.columns = ['query_ID', 'subject_ID', 'length', 'perc_id', 'q_start', 'q_end', 's_start',
                                                   's_end', 'evalue', 'q_seq', 's_seq']

                        sstart_pos = piece2_blast_df.loc[0, 's_start']

                        ### Extract 300nt harbors
                        ### Harbor from chimeric read
                        upstream_pos_int = sstart_pos
                        ###downstream_LG_primary = sequence[upstream_pos_int - 1:(sstart_pos + 300)- 1].upper()

                        downstream_LG_primary = extract_sequence(int(upstream_pos_int - 1), int((sstart_pos + 300)- 1))
                        downstream_LG_primary = downstream_LG_primary.upper()

                        ### Harbor from continous read
                        start_pos_term_RAW_LP = ((start_pos_continous_RAW - 150) + cigar_term - 1)
                        upstream_pos_term = ((start_pos_continous_RAW - 150) + cigar_term) - 299
                        ##upstream_LG_primary = sequence[upstream_pos_term - 1:start_pos_term_RAW_LP].upper()
                        upstream_LG_primary = extract_sequence(int(upstream_pos_term - 1), int(start_pos_term_RAW_LP))
                        upstream_LG_primary = upstream_LG_primary.upper()

                        #### Purine code harbors
                        upstream_LG = primary2purine(str(upstream_LG_primary))
                        downstream_LG = primary2purine(str(downstream_LG_primary))

                        #### Blastn between harbors
                        with tempfile.NamedTemporaryFile(prefix='ext_frag_', suffix='.fa', delete=True) as upstream_LG_f, \
                                tempfile.NamedTemporaryFile(prefix='int_frag_', suffix='.fa', delete=True) as downstream_LG_f:

                            with open(upstream_LG_f.name, 'w') as f:
                                f.write(f">upstream_LG\n{upstream_LG}\n")
                            with open(downstream_LG_f.name, 'w') as f2:
                                f2.write(f">downstream_LG\n{downstream_LG}\n")

                            command = 'blastn -task blastn-short -ungapped -query {} -subject {} -qcov_hsp_perc {} -perc_identity {} -outfmt "6 qseqid sseqid length pident qstart qend sstart send evalue qseq sseq"'.format(
                                upstream_LG_f.name, downstream_LG_f.name, args.cov, float(args.perc))
                            res_blast = subprocess.run(command, shell=True, capture_output=True, text=True)
                            output_string = res_blast.stdout

                            if output_string.strip():
                                output_file_like_object = StringIO(output_string)
                                #### Read the file-like object using pandas
                                result_blast_df = pd.read_csv(output_file_like_object, sep='\t', header=None)
                                result_blast_df.columns = ['query_ID', 'subject_ID', 'length', 'perc_id', 'q_start', 'q_end', 's_start',
                                                           's_end', 'evalue', 'q_seq', 's_seq']
                                ###print(result_blast_df)
                                if strand == "same":
                                    ### Keep only +/+ and -/- alignments:
                                    specific_orientation = result_blast_df[((result_blast_df['q_start'] < result_blast_df['q_end']) & (result_blast_df['s_start'] < result_blast_df['s_end'])) | \
                                                                           ((result_blast_df['q_start'] > result_blast_df['q_end']) & (result_blast_df['s_start'] > result_blast_df[ 's_end']))]
                                elif strand == "opposite":
                                    ### Keep only +/- and -/+ alignments:
                                    specific_orientation = result_blast_df[((result_blast_df['q_start'] < result_blast_df['q_end']) & (result_blast_df['s_start'] > result_blast_df['s_end'])) | \
                                                                           ((result_blast_df['q_start'] > result_blast_df['q_end']) & (result_blast_df['s_start'] < result_blast_df['s_end']))]
                                specific_orientation = result_blast_df

                                if not specific_orientation.empty:
                                    #specific_orientation = specific_orientation.iloc[specific_orientation['length'].idxmax()].to_frame().T
                                    specific_orientation = specific_orientation.loc[
                                        specific_orientation['length'].idxmax()].to_frame().T
                                    appended_df = None
                                    df_result2 = pd.DataFrame()
                                    for row_tuple in specific_orientation.itertuples(index=False):
                                        ###upstream_pos_int = cont_read - 299
                                        ##upstream_pos_int = upstream_pos_int
                                        ###upstream_pos_term = terminal_end - 299
                                        ###upstream_pos_term = upstream_pos_term

                                        end_pos_internal = (sstart_pos + 300)- 1
                                        end_pos_terminal = start_pos_term_RAW_LP

                                        satisfied_homologs = calc_ligdist(row_tuple, upstream_pos_term, end_pos_terminal,
                                                                          upstream_pos_term, end_pos_internal, appended_df)

                                        if satisfied_homologs is not None or not satisfied_homologs.empty:
                                            first_row_values = satisfied_homologs.loc[
                                                0, ['end_harbor', 'start_harbor2']]

                                            #### Concatenating the values with an underscore
                                            out_file = '_'.join(first_row_values.astype(str))

                                            new_column = pd.DataFrame(
                                                {'LigP_ID': [f"LigP_{out_file}"] * len(satisfied_homologs)})

                                            #### Concatenate the new DataFrame with the original DataFrame, placing the new column at the beginning
                                            satisfied_homologs = pd.concat([new_column, satisfied_homologs], axis=1)
                                            print(f"SATIFIED\n{satisfied_homologs}")
                                            # print(satisfied_homologs)
                                            #### Save as dataframe
                                            satisfied_homologs.to_csv(
                                                f"{out_path}dotplots_opposite/proximal_LigP_finder_output.tsv",
                                                mode='a',
                                                header=False, index=False, sep='\t')

##chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
chromosomes = ['chr2']

for chromosome in chromosomes:
    folder_chr = f"{out_folder}{chromosome}"
    if check_file(f"{folder_chr}/dotplots_opposite/proximal_LigP_finder_output.tsv") == True:
        print(f"Results for {chromosome} have been found! Skipping...")
        pass
    else:
        print(f"Running analysis for {chromosome}...")
        if check_file(f"{folder_chr}/dotplots_opposite/proximal_LigP_finder_output.tsv") == False:
            print(f"Homology graphs opposite direction torwards LP - proximal_spotchecking.py")
            # if os.path.exists(f"{folder_chr}/dotplots_opposite"):
            #     os.rmdir(f"{folder_chr}/dotplots_opposite")
            # os.mkdir(f"{folder_chr}/dotplots_opposite")
            if check_file(f"{folder_chr}/all_chim_reads_R2.fastq") == False:
                print("Creating all_chim_reads fastq file...")
                subprocess.run('cut -f1 {}/chimeric_reads_rmDUP.tsv > {}/listIDs.lst'.format(folder_chr, folder_chr), shell = True)
                subprocess.run('seqtk subseq {} {}/listIDs.lst > {}/all_chim_reads_R1.fastq'.format(str(args.read1), folder_chr, folder_chr, folder_chr), shell = True)
                subprocess.run('seqtk subseq {} {}/listIDs.lst > {}/all_chim_reads_R2.fastq'.format(str(args.read2), folder_chr, folder_chr, folder_chr), shell = True)
            import proximal_spotchecking
            proximal_spotchecking.process_prox(folder_chr, f"{args.genome}", chromosome)
            print("Done")
print(f"LP_finder has been done successfully!")