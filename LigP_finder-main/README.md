# LigP_finder_v2
LigP finder has been developed to detect ligation points (LPs) using micro-C data. For all LPs detected, the code searched for purine homology between harbors, which are defined by 300nt distal sequences nearby the predicted LP.

## Requirements

```
python 3.12

Packages: 
biopython	1.83
et-xmlfile	1.1.0
numpy	1.26.4
openpyxl	3.1.2
pandas	2.2.2
pip	23.2.1
pysam	0.22.0
python-dateutil	2.9.0.post0
pytz	2024.1
six	1.16.0
tqdm	4.66.2
tzdata	2024.1

Softwares:
bowtie2 2.5.2
samtools 1.19.2
bwa 0.7.17

```

The algorithm has been tested with public micro-C data deposited on SRA database from NCBI, under accession number SRR12625672 https://www.ncbi.nlm.nih.gov/sra/?term=SRR12625672

The data corresponds to 42.95M paired-end reads with 150nt each. The data has high quality and corresponds to HUDEP cell line.

In LigP finder v2, it performs a whole-genome alignment. So the user must provide the fasta file with the chromosomes' sequences with `--genome`parameter. The preliminary results were obtained with hg38, unmasked version: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

Download it into the `data/` folder, as well as fastq reads that will use (we used SRR12625672 and SRR12625674).

Alongside with the genome, you must use the parameters `--read1`and `--read2` with the paired-end reads from micro-C data.

By default, the folder LP_finder_output will be created with the output. In this folder, you'll have 25 folder, one per chromosome. Finally, in each chromosome's folder you will have intermediate files, such as "chimeric_rads.tsv", raw sam file, and others. Your result will be found in the folders dotplots_same and dotplots_opposite (for each chromosome!).


You can see the usage options of LigP finder by running `python LP_finder_whole_genome.py --help`

```
This script has been designed to identify ligation points, extract harbors, and find homologs between them.

options:
  -h, --help        show this help message and exit

Required arguments:
  --genome          Genome in fasta
  --read1           Fastq R1 from paired-end
  --read2           Fastq R2 from paired-end

Optional arguments:
  --harbor_length   Harbor length in nucleotides. Default = 300
  --cov             Minimum percentage length of the harbor to be an homolog. Default = 10 (10% of 300, 30nt per 300nt harbor)
  --perc            Minimum homology percentage between homologs. Default = 80
  --out             Folder name to make the output files. Default = LP_finder_output
  --dotplot         Activate the creation of 50 dotplots per orientation
```

IMPORTANT: 
2. Intermediate files can be heavy (sam and bam files), be aware to have at least 100Gb available




## Proximal and distal controls

- Proximal makes an artificial chimeric from the continous read, following how its counterpart (real chimeric read) has been splitted into two pieces. 
- Distal makes an artifical chimeric read mixing one piece from the chimeric read with another piece from the continous read.

Proximal and distal control work pretty much as the "whole_genome" script. So all the steps since the whole-genome alignment with bowtie2 to the identification of LPs per chromosome are not performed again. These two control codes will simply use the LPs that were identified in each chromosome and make harbors in a different way.

In other words, the two controls will search for homologs using different harbors sequences, that's why we don't need to rerun all the steps before LP identificaton (table chimeric_reads.tsv).

Just be aware of running `LP_finder_whole_genome.py` before the controls to have all necessary intermediate files.
