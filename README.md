# reimagined-octo-funicular
Requirements
-Python
-Biopython
-HMMER3
-MEME-Suite
-Pfam-A database

Instructions for running lanthipeptide genome mining pipeline:

Note: replace text between < > with what is appropriate for you

Preparing genbank files
1) Download genbank and nucleotide fasta files from NCBI
2) Move file prepare_genbank/prepare-genbank-files.py to directory with the genbank and nucleotide fasta files
3) Edit the range in prepare-genbank-files.py to reflect the number of files you are preparing
4) Execute prepare-genbank-files.py
5) Move file prepare_genbank/index_gb.py to the directory with the genbank files
6) Execute index_gb.py

Extracting potential lanthipeptide biosynthetic gene clusters
1) Edit get_lanC.py so fh = open(<list-of-lanC-like-protein-accession-numbers.txt>) and fh2 = open(<directory_with_genbank_files>/ref_seq-acc.csv')
2) Execute get_lanC.py
3) Edit lanthipeptide_gb_mining.py so variables at the beginning of file match your computer
4) Execute lanthipeptide_gb_mining.py
5) Execute split_classes.py (may need to alter Pfam family IDs if versions are different in your Pfam database)
6) Edit pull-classes.py so fh1 = open(<class of lanthipeptide you are analyzing>) and fh2 = open('lanC_like-acc.csv')
7) Execute pull-classes.py, note results are put to stdout
8) Edit pull-classes.py so fh2 = open('lanC-like-peps.csv')
9) Execute pull-classes.py, note results are put to stdout

Identifying potential precursor peptides in class I clusters
1) Edit make-fasta.py so fh = open(<classI-potential precursor peptides file.csv from previous step>)
2) Execute make-fasta.py, note results are put to stdout
3) Execute fimo classI_leader-meme.txt <classI potential precursor peptides fasta file>
4) Copy fimo.tsv from created directory to current directory
5) Execute get_params_classI.py <classI-potential precursor peptides file.csv> <classI-peps-params.csv>
6) Edit select-features.py so fh2 = open(<classI-peps-params.csv>)
7) Execute select-features.py
8) Execute classI_svm_classify.py
9) Execute fimo-compile.py
10) Open classI-peps-params.csv with excel
11) Delete the columns with the amino acid pairs feature columns
12) Open classI-peps-fimo-scores.csv with excel
13) Copy fimo scores and paste into the params file
14) Open classI-peps-classification.csv 
15) Copy classifications and paste into params file
16) Calculate score for peptides (+5 for svm classification, +5 for match with leader meme, +2 for core pI less than 9, +2 for 2 or more Cys in core, +1 if leader ends in GG)
17) If peptide score is greater than 10, it is predicted to be a precursor peptide
  


