import subprocess
import os
import sys
import re
import tempfile
import argparse
from io import StringIO
from Bio import SeqIO
import pandas as pd
pd.set_option('display.max_columns', None)

sys.path.insert(1, 'scripts/')

parser = argparse.ArgumentParser(description='This script has been designed to identify ligation points, extract harbors, and find homologs between them.')

required = parser.add_argument_group('Required arguments')
required.add_argument('--genome', help='Genome in fasta', required=True, type=str, metavar="")
required.add_argument('--read1', help='Fastq R1 from paired-end', required=True, type=str, metavar="")
required.add_argument('--read2', help='Fastq R2 from paired-end', required=True, type=str, metavar="")

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('--harbor_length', help='Harbor length in nucleotides. Default = 300', required=False, type=int, default=300, metavar="")
optional.add_argument('--cov', help='Minimum percentage length of the harbor to the homologs. Default = 10', required=False, type=float, default=10, metavar="")
optional.add_argument('--perc', help='Minimum homology percentage between homologs. Default = 80', required=False, type=int, default=80, metavar="")
optional.add_argument('--out', help='Folder name to make the output files. Default = LP_finder_output', required=False, type=str, default="LP_finder_output", metavar="")
optional.add_argument('--dotplot', help='Activate the creation of 50 dotplots per orientation', required=False, action = 'store_true')

args = parser.parse_args()

##print(args.dotplot)
def check_file(f):
    if os.path.exists(f):
        if os.stat(str(f)).st_size > 0:
            return True
        else:
            return False
    else:
        return False

genome_file = args.genome
out_folder = str(f"{args.out}/")

if check_file(f"{args.out}") == False:
    os.mkdir(f"{args.out}")

def complementar(seq):
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
    #     if out_df is None:
    #         out_df = result_blast
    #     else:
    #         out_df = pd.concat([out_df, out_df], axis=1)
    #         print(out_df)
    # if out_df is not None:
        return result_blast

def primary2purine(seq):
    purine = re.sub(r'a|A', 'G', str(seq))
    purine = re.sub(r't|T', 'C', str(purine))
    return purine

def homologs_blast(row, sequence, strand, load_pbar):
    load_pbar.update(1)
    start_pos_int_RAW = row[3]
    # if start_pos_int == 2840597:
    start_pos_term_RAW = row[6]
    #print(start_pos_term_RAW)
    #### Getting read sequence
    # chim_reads = pd.read_csv(f"data/chimeric_reads.tsv", sep='\t')

    # #### Get read split into the two positions
    # read_ID = chim_reads[(chim_reads['internal_arm'] == start_pos_int_RAW) & (chim_reads['terminal_arm'] == start_pos_term_RAW)].iloc[0]
    # read_ID = read_ID['read_id']

    ####Get raw read sequence
    # read_sequence = extract_sequence_by_id("data/SGreads_merged.fastq", read_ID)
    # read_sequence_rev = complementar(read_sequence)
    ### Internal arm length
    cigar_int = row[4]
    cigar_int = re.search(r'(\d+)M', cigar_int)
    cigar_int = int(cigar_int.group(1))

    ### Terminal arm length
    cigar_term = row[7]
    cigar_term = re.search(r'(\d+)M', cigar_term)
    cigar_term = int(cigar_term.group(1))

    if start_pos_int_RAW > start_pos_term_RAW:
        start_pos_int_RAW, start_pos_term_RAW = start_pos_term_RAW, start_pos_int_RAW
        cigar_int, cigar_term = cigar_term, cigar_int

    # ### Get read pieces
    # read_internal = read_sequence[:cigar_int]
    # read_internal_rev = complementar(read_internal)
    # # print(f"Piece 1:\n{read_internal}")
    #
    # read_terminal = read_sequence[cigar_int:]
    # read_terminal_rev = complementar(read_terminal)



    # ### 30nt flanking region LP1
    start_pos_int_ad = start_pos_int_RAW + cigar_int - 1
    # upstream_LP1 = start_pos_int_ad - 29
    # downstream_LP1 = start_pos_int_ad + 29
    # up_seq_LP1 = sequence[upstream_LP1 - 1:start_pos_int_ad].upper()
    # down_seq_LP1 = sequence[start_pos_int_ad - 1:downstream_LP1].upper()
    # up_seq_LP1_rev = complementar(up_seq_LP1)
    # downstream_LP1_rev = complementar(down_seq_LP1)

    # # ### 30nt flanking region LP2
    start_pos_term_ad = start_pos_term_RAW + cigar_term - 1
    # upstream_LP2 = start_pos_term_ad - 29
    # up_seq_LP2 = sequence[upstream_LP2 - 1:start_pos_term_ad].upper()
    # downstream_LP2 = start_pos_term_ad + 29
    # down_seq_LP2 = sequence[start_pos_term_ad - 1:downstream_LP2].upper()
    # ### Get reverse strand for flanking regions
    # up_seq_LP2_rev = complementar(up_seq_LP2)
    # downstream_LP2_rev = complementar(down_seq_LP2)
    #
    # up_comp_rev_LP1 = complementar_rev(up_seq_LP1_rev)
    # up_comp_rev_LP2 = complementar_rev(up_seq_LP2_rev)



    ### Extract 300nt harbors
    upstream_pos_int = start_pos_int_ad - 299
    upstream_LG_primary = sequence[upstream_pos_int - 1:start_pos_int_ad].upper()

    upstream_pos_term = start_pos_term_ad - 299
    downstream_LG_primary = sequence[upstream_pos_term - 1:start_pos_term_ad].upper()

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

        if res_blast.returncode == 0:
            #### Convert the output string into a file-like object
            output_string = res_blast.stdout

            output_file_like_object = StringIO(output_string)
            #### Read the file-like object using pandas
            result_blast_df = pd.read_csv(output_file_like_object, sep='\t', header=None)
            if not result_blast_df.empty:
                result_blast_df.columns = ['query_ID', 'subject_ID', 'length', 'perc_id', 'q_start', 'q_end', 's_start',
                                           's_end', 'evalue', 'q_seq', 's_seq']

                if strand == "same":
                    ### Keep only +/+ and -/- alignments:
                    specific_orientation = result_blast_df[((result_blast_df['q_start'] < result_blast_df['q_end']) & (result_blast_df['s_start'] < result_blast_df['s_end'])) | \
                                                           ((result_blast_df['q_start'] > result_blast_df['q_end']) & (result_blast_df['s_start'] > result_blast_df[ 's_end']))]
                elif strand == "opposite":
                    ### Keep only +/- and -/+ alignments:
                    specific_orientation = result_blast_df[((result_blast_df['q_start'] < result_blast_df['q_end']) & (result_blast_df['s_start'] > result_blast_df['s_end'])) | \
                                                           ((result_blast_df['q_start'] > result_blast_df['q_end']) & (result_blast_df['s_start'] < result_blast_df['s_end']))]
                if not specific_orientation.empty:
                    specific_orientation = specific_orientation.iloc[specific_orientation['length'].idxmax()].to_frame().T
                    appended_df = None
                    df_result2 = pd.DataFrame()
                    for row_tuple in specific_orientation.itertuples(index=False):
                        # harbor1_start = row_tuple[4]
                        # harbor1_end = row_tuple[5]
                        # harbor2_start = row_tuple[6]
                        # harbor2_end = row_tuple[7]
                        #
                        # harbor1_marked = insert_brackets(upstream_LG_primary, harbor1_start, harbor1_end, strand, False)
                        # qseq = upstream_LG_primary[harbor1_start:harbor1_end - 1]
                        # harbor1_marked_rev = complementar(harbor1_marked)


                        # if strand == "opposite":
                        #     sseq = downstream_LG_primary[harbor2_end: harbor2_start - 1]
                        #     ### ADD COMP REV TOO
                        #     harbor2_marked = insert_brackets(downstream_LG_primary, harbor2_start, harbor2_end, strand, True)
                        # else:
                        #     sseq = downstream_LG_primary[harbor2_start:harbor2_end - 1]
                        #     harbor2_marked = insert_brackets(downstream_LG_primary, harbor2_start, harbor2_end, strand,False)
                        # harbor2_marked_rev = complementar(harbor2_marked)
                        # harbor2_marked_rev_comp = complementar_rev(harbor2_marked_rev)
                        # harbor2_marked_rev_comp = re.sub(r'>>|<<', '', str(harbor2_marked_rev_comp))

                        # if check_file(f"{strand}_homolog_primary.txt") == False:
                        #     with open(f"{strand}_homolog_primary.txt", 'w') as f_out:
                        #         if strand == "same":
                        #             ###PRIMARY CODE
                        #             f_out.write(f"Same orientation results\n\n")
                        #         elif strand == "opposite":
                        #             f_out.write(f"opposite orientation results\n\n")
                        #         ### Raw read
                        #         f_out.write(f"Raw read sequence: {read_ID}\n5'{read_sequence}3'\n3'{read_sequence_rev}5'\n\n")
                        #         ### Read with LP marked
                        #         f_out.write(f"Read with LP marked:\n5'{read_internal}O{read_terminal}3'\n3'{read_internal_rev}O{read_terminal_rev}5'\n\n\n\n\n\n")
                        #
                        #         ### Piece 1
                        #         f_out.write(f"Piece 1\n5'{read_internal}O3'\n3'{read_internal_rev}O5'\n\n")
                        #         # ### Flanking 30nt
                        #         f_out.write(f"flanking LP1\n5'{up_seq_LP1}O{down_seq_LP1}3'\n3'{up_seq_LP1_rev}O{downstream_LP1_rev}5'\n\n")
                        #         ### Harbors
                        #         f_out.write(f"Harbor1_homolog_{upstream_pos_int}_{start_pos_int_ad}\n5'{harbor1_marked}O3'\n3'{harbor1_marked_rev}O5'\n\n\n\n\n\n")
                        #
                        #         ### Piece 2
                        #         f_out.write(f"Piece 2\n5'O{read_terminal}3'\n3'O{read_terminal_rev}5'\n\n")
                        #         ### Flanking 30nt
                        #         f_out.write(f"flanking LP2\n5'{up_seq_LP2}O{down_seq_LP2}3'\n3'{up_seq_LP2_rev}O{downstream_LP2_rev}5'\n<-{'--' * int(27)}\nc{up_comp_rev_LP2}\n")
                        #         ### Harbors
                        #         f_out.write(f"\nHarbor2_homolog_{upstream_pos_term}_{start_pos_term_ad}\n5'{harbor2_marked}O3'\n3'{harbor2_marked_rev}O5'\n<{'-' * int(270)}\nc{harbor2_marked_rev_comp}\n\n\n")
                        #
                        #         if strand == "opposite":
                        #             f_out.write(
                        #                 f"Homologs alignment - purine\n\n5'{qseq}3'\n5'{complementar_and_rev(sseq)}3'")
                        #         else:
                        #             f_out.write(
                        #                 f"Homologs alignment - purine\n\n5'{qseq}3'\n5'{sseq}3'")
                        #
                        #
                        #
                        #
                        #         ###PURINE CODE
                        #         f_out.write(f"\n\n\n\n Purine code \n\n")
                        #         ### Raw read
                        #         f_out.write(f"Raw read sequence: {primary2purine(read_ID)}\n5'{primary2purine(read_sequence)}3'\n3'{primary2purine(read_sequence_rev)}5'\n\n")
                        #         ### Read with LP marked
                        #         f_out.write(f"Read with LP marked:\n5'{primary2purine(read_internal)}O{primary2purine(read_terminal)}3'\n3'{primary2purine(read_internal_rev)}O{primary2purine(read_terminal_rev)}5'\n\n\n\n\n\n")
                        #
                        #         ### Piece 1
                        #         f_out.write(f"Piece 1\n5'{primary2purine(read_internal)}O3'\n3'{primary2purine(read_internal_rev)}O5'\n\n")
                        #         # ### Flanking 30nt
                        #         f_out.write(f"flanking LP1\n5'{primary2purine(up_seq_LP1)}O{primary2purine(down_seq_LP1)}3'\n3'{primary2purine(up_seq_LP1_rev)}O{primary2purine(downstream_LP1_rev)}5'\n\n")
                        #         ### Harbors
                        #         f_out.write(f"Harbor1_homolog_{primary2purine(upstream_pos_int)}_{primary2purine(start_pos_int_ad)}\n5'{primary2purine(harbor1_marked)}O3'\n3'{primary2purine(harbor1_marked_rev)}O5'\n\n\n\n\n\n")
                        #
                        #         ### Piece 2
                        #         f_out.write(f"Piece 2\n5'O{primary2purine(read_terminal)}3'\n3'O{primary2purine(read_terminal_rev)}5'\n\n")
                        #         ### Flanking 30nt
                        #         f_out.write(f"flanking LP2\n5'{primary2purine(up_seq_LP2)}O{primary2purine(down_seq_LP2)}3'\n3'{primary2purine(up_seq_LP2_rev)}O{primary2purine(downstream_LP2_rev)}5'\n<-{'--' * int(27)}\nc{primary2purine(up_comp_rev_LP2)}\n")
                        #         ### Harbors
                        #         f_out.write(f"\nHarbor2_homolog_{primary2purine(upstream_pos_term)}_{primary2purine(start_pos_term_ad)}\n5'{primary2purine(harbor2_marked)}O3'\n3'{primary2purine(harbor2_marked_rev)}O5'\n<{'-' * int(270)}\nc{primary2purine(harbor2_marked_rev_comp)}\n\n")
                        #
                        #         if strand == "opposite":
                        #             f_out.write(f"Homologs alignment - purine\n\n5'{primary2purine(qseq)}3'\n5'{primary2purine(complementar_and_rev(sseq))}3'")
                        #         else:
                        #             f_out.write(f"Homologs alignment - purine\n\n5'{primary2purine(qseq)}3'\n5'{primary2purine(sseq)}3'")
                        upstream_pos_int = start_pos_int_ad - 299
                        upstream_pos_term = start_pos_term_ad - 299
                        satisfied_homologs = calc_ligdist(row_tuple, upstream_pos_int, start_pos_int_ad,
                                                          upstream_pos_term, start_pos_term_ad, appended_df)

                        if satisfied_homologs is not None or not satisfied_homologs.empty:

                            first_row_values = satisfied_homologs.loc[0, ['end_harbor', 'start_harbor2']]

                            #### Concatenating the values with an underscore
                            out_file = '_'.join(first_row_values.astype(str))

                            new_column = pd.DataFrame({'LigP_ID': [f"LigP_{out_file}"] * len(satisfied_homologs)})

                            #### Concatenate the new DataFrame with the original DataFrame, placing the new column at the beginning
                            satisfied_homologs = pd.concat([new_column, satisfied_homologs], axis=1)

                            #### Save as dataframe
                            satisfied_homologs.to_csv(f"{folder_chr}/dotplots_{strand}/LigP_finder_output_{strand}.tsv", mode='a',
                                                      header=False, index=False, sep='\t')

                    png_files = [file for file in os.listdir(f"{folder_chr}/dotplots_{strand}") if file.endswith('.png')]
                    num_dotplots = len(png_files)
                    if num_dotplots <= 50:
                        ######## Dense dotplots
                        command = 'blastn -task blastn-short -query {} -subject {} -qcov_hsp_perc 2.5 -perc_identity 95 -ungapped -outfmt "6 qseqid sseqid length pident qstart qend sstart send evalue qseq sseq"'.format(
                            upstream_LG_f.name, downstream_LG_f.name)
                        res_blast = subprocess.run(command, shell=True, capture_output=True, text=True)
                        if res_blast.returncode == 0:
                            # Convert the output string into a file-like object
                            output_string = res_blast.stdout
                            output_file_like_object = StringIO(output_string)
                            # Read the file-like object using pandas
                            result_blast_df = pd.read_csv(output_file_like_object, sep='\t', header=None)
                            result_blast_df.columns = ['query_ID', 'subject_ID', 'length', 'perc_id', 'q_start',
                                                       'q_end',
                                                       's_start',
                                                       's_end', 'evalue', 'q_seq', 's_seq']

                            if strand == "same":
                                ### Keep only +/+ and -/- alignments:
                                specific_orientation = result_blast_df[((result_blast_df['q_start'] < result_blast_df[
                                    'q_end']) & (result_blast_df['s_start'] < result_blast_df['s_end'])) | (
                                                                                   (result_blast_df[
                                                                                        'q_start'] >
                                                                                    result_blast_df[
                                                                                        'q_end']) & (
                                                                                           result_blast_df[
                                                                                               's_start'] >
                                                                                           result_blast_df[
                                                                                               's_end']))]
                            elif strand == "opposite":
                                ### Keep only +/- and -/+ alignments:
                                specific_orientation = result_blast_df[((result_blast_df['q_start'] < result_blast_df[
                                    'q_end']) & (result_blast_df['s_start'] > result_blast_df['s_end'])) | (
                                                                                   (result_blast_df[
                                                                                        'q_start'] <
                                                                                    result_blast_df[
                                                                                        'q_end']) & (
                                                                                           result_blast_df[
                                                                                               's_start'] >
                                                                                           result_blast_df[
                                                                                               's_end']))]

                            appended_df = None
                            df_result2 = pd.DataFrame()

                            for row_tuple_dense in specific_orientation.itertuples(index=False):
                                satisfied_homologs = calc_ligdist(row_tuple_dense, upstream_pos_int, start_pos_int_ad,
                                                                  upstream_pos_term, start_pos_term_ad, appended_df)

                                if satisfied_homologs is not None or not satisfied_homologs.empty:
                                    first_row_values = satisfied_homologs.loc[0, ['end_harbor', 'start_harbor2']]
                                    #### Concatenating the values with an underscore
                                    out_file = '_'.join(first_row_values.astype(str))

                                    new_column = pd.DataFrame(
                                        {'LigP_ID': [f"LigP_{out_file}"] * len(satisfied_homologs)})

                                    #### Concatenate the new DataFrame with the original DataFrame, placing the new column at the beginning
                                    satisfied_homologs = pd.concat([new_column, satisfied_homologs], axis=1)

                                    concat_columns = pd.concat(
                                        [satisfied_homologs.iloc[:, 0], satisfied_homologs.iloc[:, 5:9]], axis=1)
                                    df_result2 = pd.concat([df_result2, concat_columns], ignore_index=True)


                            with tempfile.NamedTemporaryFile(prefix='blast_out', suffix='.tsv',
                                                             delete=True) as blast_table:
                                df_result2.to_csv(f"{blast_table.name}", header=True, index=False, sep='\t')
                                subprocess.run(["Rscript", "--vanilla", "scripts/dotplot.R", blast_table.name, out_file,
                                                f"{folder_chr}/dotplots_{strand}"])


####Create bowtie2 index
if check_file(str(f"{out_folder}hg38_unmasked.1.bt2")) == False:
    print(f"Creating bowtie2 index...")
    subprocess.call(['bowtie2-build', f"{args.genome}", f"{out_folder}hg38_unmasked", '--threads', '8'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print(f"Done")

####Perform end to end alignment 100% homology
if check_file(str(f"{out_folder}microC.sam")) == False and check_file(str(f"{out_folder}microC_SGreads.sam")) == False:
    print(f"Performing bowtie2 alignment...")
    subprocess.call(['bowtie2', '-x', f"{out_folder}hg38_unmasked", '-1', f"{args.read1}", '-2', f"{args.read2}", '-S', f"{out_folder}microC.sam", '--threads', '8'])
    print(f"Done")

print("Performing analysis chromosome per chromosome...")
##chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
chromosomes = ["chr2"]
### Define input SAM file
input_sam_file = f"{out_folder}microC.sam"

### Iterate through each chromosome
for chromosome in chromosomes:
    folder_chr = f"{out_folder}{chromosome}"
    if check_file(f"{folder_chr}/dotplots_opposite/LigP_finder_output_same.tsv") == True:
        print(f"Results for {chromosome} have been found! Skipping...")
        pass
    else:
        print(f"Running analysis for {chromosome}...")
    if check_file(f"{folder_chr}") == False:
        os.mkdir(f"{folder_chr}")

    chr_sam = f"{folder_chr}/{chromosome}_subset.sam"

    if check_file(f"{chr_sam}") == False:
        print(f"Creating subset from sam file for {chromosome}...")
        with open(input_sam_file, 'r') as input_file, open(chr_sam, 'w') as output_file:
            # Iterate through each line in the input SAM file
            for line in input_file:
                # Skip header lines
                if line.startswith('@'):
                    output_file.write(line)
                    continue

                # Split the SAM line into fields
                fields = line.split('\t')

                # Check if the read aligns to the current chromosome
                if fields[2] == chromosome:
                    # Write the line to the output SAM file
                    output_file.write(line)
        print("Done!")

    ####Extract singleton reads
    if check_file(str(f"{folder_chr}/{chromosome}_SGreads.sam")) == False:
        print(f"Extract singleton reads")
        subprocess.call(['samtools', 'view', '-@', '8', '-f', '4', '-F', '8', f"{chr_sam}", '--output-fmt', 'SAM', '-o', f"{folder_chr}/{chromosome}_SGreads.sam"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"Done")

    ####Take only R1 singletons
    if check_file(str(f"{folder_chr}/{chromosome}_SGreads_R1.sam")) == False:
        print(f"Take only R1 singletons")
        subprocess.run('head -26 {} | grep "@" > {}/temp.txt'.format(chr_sam, folder_chr), shell = True)
        subprocess.run('cat {}/temp.txt {}/{}_SGreads.sam > {}/tmp.sam && mv {}/tmp.sam {}/{}_SGreads.sam'.format(folder_chr,folder_chr,chromosome,
                                                                                                                  folder_chr,
                                                                                                                  folder_chr,
                                                                                                                  folder_chr,
                                                                                                                  chromosome), shell = True)

        subprocess.run('samtools view -@ 8 -f 64 {}/{}_SGreads.sam -o {}/{}_SGreads_R1.sam'.format(folder_chr, chromosome,
                                                                                                     folder_chr, chromosome), shell = True)

    ####Take only R2 singletons
    if check_file(str(f"{folder_chr}/{chromosome}_SGreads_R2.sam")) == False:
        print(f"Take only R2 singletons")
        subprocess.run('samtools view -@ 8 -f 128 {}/{}_SGreads.sam -o {}/{}_SGreads_R2.sam'.format(folder_chr, chromosome,
                                                                                                     folder_chr, chromosome), shell = True)
        print(f"Done")

    ####Extract singleton reads R1 and R2 merged to a fastq file
    if check_file(str(f"{folder_chr}/{chromosome}_SGreads_R1.fastq")) == False:
        print("Extract singleton R1 reads and create fastq file")
        subprocess.run('cat {}/temp.txt {}/{}_SGreads_R1.sam > {}/tmp2.txt && mv {}/tmp2.txt {}/{}_SGreads_R1.sam'.format(folder_chr, folder_chr, chromosome,
                                                                                                                              folder_chr, folder_chr, folder_chr, chromosome), shell = True)
        subprocess.run('samtools view -@ 8 -bS {}/{}_SGreads_R1.sam -o {}/{}_SGreads_R1.bam'.format(folder_chr, chromosome,
                                                                                                    folder_chr, chromosome), shell = True)
        subprocess.run('samtools sort -@ 8 -o {}/{}_SGreads_R1_ST.bam {}/{}_SGreads_R1.bam'.format(folder_chr, chromosome,
                                                                                                    folder_chr, chromosome), shell = True)
        subprocess.run('samtools fastq {}/{}_SGreads_R1_ST.bam -o {}/{}_SGreads_R1.fastq'.format(folder_chr, chromosome,
                                                                                                    folder_chr, chromosome), shell = True)

    if check_file(str(f"{folder_chr}/{chromosome}_SGreads_R2.fastq")) == False:
        print("Extract singleton R2 reads and create fastq file")
        subprocess.run('cat {}/temp.txt {}/{}_SGreads_R2.sam > {}/tmp2.txt && mv {}/tmp2.txt {}/{}_SGreads_R2.sam'.format(folder_chr, folder_chr, chromosome,
                                                                                                                          folder_chr, folder_chr, folder_chr, chromosome), shell=True)
        subprocess.run('samtools view -@ 8 -bS {}/{}_SGreads_R2.sam -o {}/{}_SGreads_R2.bam'.format(folder_chr, chromosome,
                                                                                         folder_chr, chromosome),shell=True)
        subprocess.run('samtools sort -@ 8 -o {}/{}_SGreads_R2_ST.bam {}/{}_SGreads_R2.bam'.format(folder_chr, chromosome,
                                                                                        folder_chr, chromosome),shell=True)
        subprocess.run('samtools fastq {}/{}_SGreads_R2_ST.bam -o {}/{}_SGreads_R2.fastq'.format(folder_chr, chromosome,
                                                                                     folder_chr, chromosome),shell=True)
        os.remove(f"{folder_chr}/temp.txt")

    if check_file(str(f"{folder_chr}/{chromosome}_SGreads_merged.fastq")) == False:
        subprocess.run('cat {}/{}_SGreads_R1.fastq {}/{}_SGreads_R2.fastq > {}/SGreads_merged.fastq'.format(folder_chr, chromosome,
                                                                                                        folder_chr, chromosome, folder_chr), shell = True)
        print(f"Done")

    ####Align fastq singleton reads as single-end file using bwa
    if check_file(str(f"{folder_chr}/{chromosome}_bwa_band.sam")) == False:
        print(f"Align fastq singleton reads as single-end file using bwa")
        ### Subset fasta with specific chromosome
        write_line = False
        with open(f"{args.genome}", 'r') as input_file, open(f"{folder_chr}/{chromosome}.fa", 'w') as output_file:
            for line in input_file:
                ### Check if the current line corresponds to a chromosome header
                if line.startswith('>'):
                    ### Check if the current chromosome is chromosome 1
                    if line.strip().split()[0][1:] == chromosome:
                        # Set the flag to True to start writing lines
                        write_line = True
                    else:
                        ### Set the flag to False to stop writing lines
                        write_line = False
                        continue
                ### Write the line to the output FASTA file if the flag is True
                if write_line:
                    output_file.write(line)

        subprocess.call(['bwa', 'index', f"{folder_chr}/{chromosome}.fa"], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        subprocess.call(['bwa', 'mem', '-k', '22', '-w', '1000000', '-a', '-t', '8', '-T', '0', f"{folder_chr}/{chromosome}.fa", f"{folder_chr}/SGreads_merged.fastq"],
                        stdout=open(f"{folder_chr}/{chromosome}_bwa_band.sam", 'w'), stderr=subprocess.DEVNULL)
        print(f"Done")

    ####Filter out reads with mismatches
    if check_file(str(f"{folder_chr}/{chromosome}_bwa_band_no-mismatches.sam")) == False:
        print(f"Filter out reads with mismatches")
        subprocess.run('samtools view {}/{}_bwa_band.sam | grep -w "NM:i:0" > {}/{}_bwa_band_no-mismatches.sam'.format(folder_chr, chromosome, folder_chr, chromosome), shell = True)
        print(f"Done")

    ###Get splitted reads into two fragments
    if check_file(str(f"{folder_chr}/splitted_IDs.lst")) == False:
        print(f"Get splitted reads into two fragments")
        command = "cut -f1 {}/{}_bwa_band_no-mismatches.sam | sort | uniq -c | awk '$1 == 2' | awk '{{print $2}}' > {}/splitted_IDs.lst".format(
            folder_chr, chromosome, folder_chr)
        subprocess.run(command, shell=True)

    ###Filtering out splitted reads with proximal distance
    if check_file(f"{folder_chr}/splitted_long_dist.tsv") == False:
        print(f"Filtering out splitted reads with proximal distance - distance_fragments.py")
        import distance_fragments
        distance_fragments.process_dist(folder_chr, chromosome)
        print(f"Done")

    #### Create sam file with singleton's mates
    if check_file(f"{folder_chr}/{chromosome}_SGmate.sam") == False:
        print(f"Create sam file with singleton's mates - samtools")
        subprocess.run('samtools view -@ 8 -f 8 -F 4 {} --output-fmt SAM -o {}/{}_SGmate.sam'.format(chr_sam, folder_chr, chromosome), shell = True)
        print(f"Done")

    ### Identify chimeric reads
    if check_file(f"{folder_chr}/chimeric_reads.tsv") == False:
        print(f"Identifying chimeric reads - mate_location.py")
        import mate_location
        mate_location.process(folder_chr, chromosome)
        print(f"Done")

    ### Arms homology
    if check_file(f"{folder_chr}/result_ligation_points_rmDUP.tsv") == False:
        chimeric_reads = pd.read_csv(f"{folder_chr}/chimeric_reads.tsv", sep='\t')
        chimeric_reads = chimeric_reads.drop_duplicates(subset=chimeric_reads.columns[1:7], keep='first')
        chimeric_reads.to_csv(f"{folder_chr}/tmp.tsv", header=False, index=False, sep="\t")
        subprocess.run('cut -f3-8 {}/tmp.tsv > {}/result_ligation_points_rmDUP.tsv'.format(folder_chr, folder_chr), shell = True)
        os.remove(f"{folder_chr}/tmp.tsv")
        # print(f"Measuring arms homology - terminal_extension_multi.py")
        # import terminal_extension_multi
        # if check_file(f"{folder_chr}/result_ligation_points.lst") == True:
        #     subprocess.run('ggrep -w -f {}/result_ligation_points.lst {}/chimeric_reads.tsv > {}/result_ligation_points.tsv'.format(folder_chr, folder_chr, folder_chr), shell = True)
        # if check_file(f"{folder_chr}/result_ligation_points_rmDUP.tsv") == False:
        #     subprocess.run('cut -f3-8 {}/result_ligation_points.tsv | sort | uniq > {}/result_ligation_points_rmDUP.tsv'.format(folder_chr, folder_chr), shell = True)
        # print(f"Done")

    #### Find homologs with opposite homologs
    if check_file(f"{folder_chr}/dotplots_opposite/LigP_finder_output_opposite.tsv") == False:
        print(f"Homology graphs opposite direction torwards LP - purine_homol_T90_1_opposite.py")
        if os.path.exists(f"{folder_chr}/dotplots_opposite"):
            os.rmdir(f"{folder_chr}/dotplots_opposite")
        os.mkdir(f"{folder_chr}/dotplots_opposite")
        import purine_homol_T90_1_opposite
        purine_homol_T90_1_opposite.process_opp(folder_chr, chromosome)

        print("Done")

    ### Find homologs with same direction homologs
    if check_file(f"{folder_chr}/dotplots_same/LigP_finder_output_same.tsv") == False:
        print(f"Homology graphs same direction torwards LP - purine_homol_T90_1_same.py")
        if os.path.exists(f"{folder_chr}/dotplots_same"):
            os.rmdir(f"{folder_chr}/dotplots_same")
        os.mkdir(f"{folder_chr}/dotplots_same")
        import purine_homol_T90_1_same
        purine_homol_T90_1_same.process_same(folder_chr, chromosome)
        print("Done")

    print("Subset SAM file for reads aligned to", chromosome, "has been created:", chr_sam)

print(f"LP_finder has been done successfully!")
