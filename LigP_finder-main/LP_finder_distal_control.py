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
optional.add_argument('--out', help='Folder name to make the output files. Default = LP_finder_output', required=False, type=str, default="LP_finder_output", metavar="")
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

def extract_sequence(start, end):
    fasta_file = f"data/hg38_chr8_unmasked.fa"
    fasta = pysam.FastaFile(fasta_file)

    sequence = fasta.fetch(str("chr8"), start - 1, end).upper()
    fasta.close()
    return sequence

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

def spotcheck_continous(row, strand, load_pbar, out_path, chr):
    load_pbar.update(1)
    read_ID = row[0]
    start_pos_continous_RAW = row[1]
    start_pos_term_RAW = row[6]

    #### Getting read sequence
    chim_reads = pd.read_csv(f"{out_path}/chimeric_reads_distance_distal.tsv", sep='\t')
    #### Get read split into the two positions
    read_ID = chim_reads[(chim_reads['continuous_read'] == start_pos_continous_RAW) & (chim_reads['terminal_arm'] == start_pos_term_RAW)].iloc[0]
    read_ID = read_ID['read_id']

    chim_read_sequence = extract_sequence_by_id(f"{out_path}/{chr}_SGreads_R1.fastq", read_ID)
    if chim_read_sequence is not None:
        cont_read_sequence = extract_sequence_by_id(f"{out_path}/all_chim_reads_R2.fastq", read_ID)
    elif chim_read_sequence is None:
        chim_read_sequence = extract_sequence_by_id(f"{out_path}/{chr}_SGreads_R2.fastq", read_ID)
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


    ### Get read pieces continous read
    cont_read_internal = cont_read_sequence[:cigar_int]
    cont_read_terminal = cont_read_sequence[cigar_int:]

    ### Get read pieces chimeric read
    chim_read_internal = chim_read_sequence[:cigar_term]
    chim_read_terminal = chim_read_sequence[cigar_term:]


    ### Test if continous read has 150nt long match
    with tempfile.NamedTemporaryFile(prefix='cont_read_', suffix='.fa', delete=True) as cont_read:
        with open(cont_read.name, 'w') as f:
            f.write(f">cont_read\n{cont_read_sequence}\n")
        command = 'blastn -task blastn-short -ungapped -query {} -subject {} -perc_identity 100 -qcov_hsp_perc 100 -outfmt "6 qseqid sseqid length pident qstart qend sstart send evalue qseq sseq"'.format(
            cont_read.name, f"{out_path}/{chr}.fa")
        cont_read_blast = subprocess.run(command, shell=True, capture_output=True, text=True)
        output_string = cont_read_blast.stdout
        if output_string.strip():
            output_file_like_object = StringIO(output_string)
            cont_read_blast_df = pd.read_csv(output_file_like_object, sep='\t', header=None)
            if len(cont_read_blast_df) == 1:
                cont_read_blast_df.columns = ['query_ID', 'subject_ID', 'length', 'perc_id', 'q_start', 'q_end',
                                              's_start', 's_end', 'evalue', 'q_seq', 's_seq']
                cont_length = cont_read_blast_df.loc[0, 'length']

                if cont_length == 150:
                    ### Position piece 2 from continous read
                    with tempfile.NamedTemporaryFile(prefix='piece2_', suffix='.fa', delete=True) as piece2_pos:
                        with open(piece2_pos.name, 'w') as f:
                            f.write(f">piece2\n{cont_read_terminal}\n")
                        command = 'blastn -task blastn-short -ungapped -query {} -subject {} -perc_identity 100 -qcov_hsp_perc 100 -outfmt "6 qseqid sseqid length pident qstart qend sstart send evalue qseq sseq"'.format(
                            piece2_pos.name, str(f"{out_path}/{chr}.fa"))
                        piece2_blast = subprocess.run(command, shell=True, capture_output=True, text=True)
                        #### Convert the output string into a file-like object
                        output_string = piece2_blast.stdout
                        if output_string.strip():
                            output_file_like_object = StringIO(output_string)
                            piece2_blast_df = pd.read_csv(output_file_like_object, sep='\t', header=None)
                            if len(piece2_blast_df) == 1:
                                piece2_blast_df.columns = ['query_ID', 'subject_ID', 'length', 'perc_id', 'q_start', 'q_end', 's_start',
                                                           's_end', 'evalue', 'q_seq', 's_seq']
                                sstart_pos = piece2_blast_df.loc[0, 's_start']

                                ### Extract 300nt harbors
                                upstream_pos_int = sstart_pos

                                downstream_LG_primary = extract_sequence((upstream_pos_int - 1), (sstart_pos + 150)- 1).upper()

                                start_pos_term_RAW_LP = (start_pos_term_RAW + cigar_int - 1)
                                upstream_pos_term = start_pos_term_RAW_LP - 299
                                upstream_LG_primary = extract_sequence((upstream_pos_term - 1), start_pos_term_RAW_LP).upper()
                                if 'N' not in upstream_LG_primary:
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
                                                specific_orientation = specific_orientation.loc[
                                                    specific_orientation['length'].idxmax()].to_frame().T
                                                appended_df = None
                                                df_result2 = pd.DataFrame()
                                                for row_tuple in specific_orientation.itertuples(index=False):
                                                    end_pos_internal = (sstart_pos + 300) - 1
                                                    end_pos_terminal = start_pos_term_RAW_LP

                                                    satisfied_homologs = calc_ligdist(row_tuple, upstream_pos_term,
                                                                                      end_pos_terminal,
                                                                                      upstream_pos_term, end_pos_internal,
                                                                                      appended_df)

                                                    if satisfied_homologs is not None or not satisfied_homologs.empty:
                                                        first_row_values = satisfied_homologs.loc[0, ['end_harbor', 'start_harbor2']]

                                                        #### Concatenating the values with an underscore
                                                        out_file = '_'.join(first_row_values.astype(str))

                                                        new_column = pd.DataFrame({'LigP_ID': [f"LigP_{out_file}"] * len(satisfied_homologs)})

                                                        #### Concatenate the new DataFrame with the original DataFrame, placing the new column at the beginning
                                                        satisfied_homologs = pd.concat([new_column, satisfied_homologs], axis=1)
                                                        #### Save as dataframe
                                                        satisfied_homologs.to_csv(f"{out_path}/dotplots_opposite/distal_LigP_finder_output.tsv",mode='a',header=False, index=False, sep='\t')


chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
input_sam_file = f"{out_folder}microC.sam"

for chromosome in chromosomes:
    folder_chr = f"{out_folder}{chromosome}"
    if check_file(f"{folder_chr}/dotplots_opposite/distal_LigP_finder_output.tsv") == True:
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
        if check_file(f"{folder_chr}/chimeric_reads_rmDUP.tsv") == False:
            chimeric_reads = pd.read_csv(f"{folder_chr}/chimeric_reads.tsv", sep="\t")
            chimeric_reads = chimeric_reads.drop_duplicates(subset=chimeric_reads.columns[1:7], keep='first')
            chimeric_reads.to_csv(f"{folder_chr}/chimeric_reads_rmDUP.tsv", header=False, index=False, sep="\t")
            # chimeric_reads = pd.read_csv(f"{folder_chr}/chimeric_reads.tsv", sep='\t')
            # chimeric_reads = chimeric_reads.drop_duplicates(subset=chimeric_reads.columns[1:7], keep='first')
            # chimeric_reads.to_csv(f"{folder_chr}/tmp.tsv", header=False, index=False, sep="\t")
            # subprocess.run('cut -f3-8 {}/tmp.tsv > {}/result_ligation_points_rmDUP.tsv'.format(folder_chr, folder_chr), shell = True)
            # os.remove(f"{folder_chr}/tmp.tsv")


        if check_file(f"{folder_chr}/dotplots_opposite/distal_LigP_finder_output.tsv") == False:
            print(f"Homology graphs opposite direction torwards LP - distal_spotchecking.py")
            if os.path.exists(f"{folder_chr}/dotplots_opposite"):
                os.rmdir(f"{folder_chr}/dotplots_opposite")
            os.mkdir(f"{folder_chr}/dotplots_opposite")
            if check_file(f"{folder_chr}/all_chim_reads_R2.fastq") == False:
                print("Creating all_chim_reads fastq file...")
                subprocess.run('cut -f1 {}/chimeric_reads_rmDUP.tsv > {}/listIDs.lst'.format(folder_chr, folder_chr), shell = True)
                subprocess.run('seqtk subseq {} {}/listIDs.lst > {}/all_chim_reads_R1.fastq'.format(str(args.read1), folder_chr, folder_chr, folder_chr), shell = True)
                subprocess.run('seqtk subseq {} {}/listIDs.lst > {}/all_chim_reads_R2.fastq'.format(str(args.read2), folder_chr, folder_chr, folder_chr), shell = True)
            import purine_homol_T90_1_same_distal
            purine_homol_T90_1_same_distal.process_same(folder_chr, chromosome)
            print("Done")

print(f"LP_finder has been done successfully!")
