import glob
from Bio import SeqIO

files = glob.glob("bacteria.*.genomic.gb")

stuff = SeqIO.index_db("ref_seq.idx", files, "genbank")

