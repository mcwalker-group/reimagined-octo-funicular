from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import re
import os
import time
import find_cluster
import pep_processing
import glob

acc = 'lanC-like.acc'                                                              #file with cluster accession numbers
pep_min = 30                                                                   #minimum precursor peptide length
pep_max = 120                                                      #maximum precursor peptide length
win = 7                                                                        #orfs up and down screen to examine
pfam_file = 'lanC_like-acc.csv'                                                 #file to save pfam hits
precursor_file = 'LanC-like-peps.csv'                                           #file to save potential precursor peptides
pfam_db = '/mnt/c/linux/pfam/Pfam-A.hmm'                                             #path to pfam database
E = '0.1'                                                                  #E value cutoff for hmmscan
gb_link_file = 'lanC-like_gb_ids.csv'                    #file linking genbank nucleotide accessions to protein accessions
gb_file = glob.glob('/mnt/d/bacteria_refseq/bacteria.*.genomic.gbff.bgz')             #big genbank file with sequences
index_file = '/mnt/d/bacteria_refseq/ref_seq.idx'       #index file for genbank file made with biopython SeqIO.index_db

def get_acc():
    p = 0
    proteins = []
    fh1 = open(acc)
    for line in fh1:
        proteins.append(line.strip())
        p += 1
    return proteins, p

def gb_acc_dict():
    print "getting nuc_prot link file"
    acc = {}
    fh = open(gb_link_file)
    for line in fh:
        acc[line.split(',')[0]] = line.strip().split(',')[1:]
    print "link file loaded\n"
    fh.close()
    return acc

def index_gb():
    print "indexing genbank file"
    gb_index = SeqIO.index_db(index_file, gb_file, "genbank")
    print "done indexing file\n"
    return gb_index

def get_nuc(prot_id, gb_acc):
    print 'finding nucleotide accession for ' + prot_id
    for k,v in gb_acc.iteritems():
        if prot_id in v:
            nuc_acc = k
            print "found nucleotide accession\n"
            break
        else:
            nuc_acc = "not there"
    if nuc_acc == "not there":
        print 'could not find nucleotide accession\n'
    return nuc_acc

protein_acc = open(pfam_file, 'w')
prot, p = get_acc()
output = open(precursor_file, 'w')

gb_accs = gb_acc_dict()
gb_index = index_gb()


for line in prot:
    t0 = time.time()
    found = get_nuc(line, gb_accs)
    if found != "not there":
        record = gb_index[found]
        p, phylum, species, forward_nuc, clust_coords, orf_coords, nuc_acc = find_cluster.find_cluster(record, line, p, win, E, pfam_db, protein_acc)
        print "finding precursors for: " + species
        f1 = pep_processing.get_peps(int(clust_coords[0]), forward_nuc, 'f', pep_max, pep_min)
        f1b = pep_processing.remove_redundancy(f1, 'f')
        f2 = pep_processing.get_peps(int(clust_coords[0])+1, forward_nuc[1:], 'f', pep_max, pep_min)
        f2b = pep_processing.remove_redundancy(f2, 'f')
        f3 = pep_processing.get_peps(int(clust_coords[0])+2, forward_nuc[2:], 'f', pep_max, pep_min)
        f3b = pep_processing.remove_redundancy(f3, 'f')

        rev_nuc = forward_nuc.reverse_complement()

        r1 = pep_processing.get_peps(int(clust_coords[1]), rev_nuc, 'r', pep_max, pep_min)
        r1b = pep_processing.remove_redundancy(r1, 'r')
        r2 = pep_processing.get_peps(int(clust_coords[1])-1, rev_nuc[1:], 'r', pep_max, pep_min)
        r2b = pep_processing.remove_redundancy(r2, 'r')
        r3 = pep_processing.get_peps(int(clust_coords[1])-2, rev_nuc[2:], 'r', pep_max, pep_min)
        r3b = pep_processing.remove_redundancy(r3, 'r')

        f = f1b + f2b + f3b
        r = r1b + r2b + r3b
        
        intergenic_f, intergenic_r = pep_processing.intergenic(f, r, orf_coords)
        for lan in intergenic_f:
            output.write(line + ',' + phylum + ',' + species + ',' + nuc_acc + ',' + str(lan[0]) + ',' + str(lan[1]) + ',' + 'f,' + lan[2] + '\n')
        for lan in intergenic_r:
            output.write(line + ',' + phylum + ',' + species + ',' + nuc_acc + ',' + str(lan[0]) + ',' + str(lan[1]) + ',' + 'r,' + lan[2] + '\n')
    
    t1 = time.time()
    print '\t\t timer: ' + '\t' + str(t1-t0) + '\n'
output.close()
protein_acc.close()
