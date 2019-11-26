# reimagined-octo-funicular
Instructions for running lanthipeptide genome mining pipeline:

Preparing genbank files
1) Download genbank and nucleotide fasta files from NCBI
2) Move file prepare_genbank/prepare-genbank-files.py to directory with the genbank and nucleotide fasta files
3) Edit the range in prepare-genbank-files.py to reflect the number of files you are preparing
4) Execute prepare-genbank-files.py
5) Move file prepare_genbank/index_gb.py to the directory with the genbank files
6) Execute index_gb.py

Extracting potential lanthipeptide biosynthetic gene clusters
1) Edit get_lanC.py so fh = open('list-of-lanC-like-protein-accession-numbers.txt') and fh2 = open('directory_with_genbank_files/ref_seq-acc.csv')
2) Execute get_lanC.py
3) Edit lanthipeptide_gb_mining.py so variables at the beginning of file match your computer
4) Execute lanthipeptide_gb_mining.py
5) Execute split_classes.py (may need to alter Pfam family IDs if versions are different in your Pfam database)
6) Edit pull_classes.py 

Identifying potential precursor peptides in class I clusters
1)
