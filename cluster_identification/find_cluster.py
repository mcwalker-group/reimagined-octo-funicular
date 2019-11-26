#finding biosynthetic cluster

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

def find_cluster(record, accession, p, win, E, pfam_db, protein_acc):
    i = 0
    for line in record.features:
        if line.type == "CDS":
            if "protein_id" in line.qualifiers:
                if line.qualifiers["protein_id"][0] == accession:
                    clust = {}
                    aa_out = open('clust.faa', 'w')
                    print "finding cluster: " + line.qualifiers["protein_id"][0] + ' (' + str(p) + ")"
                    p = p - 1
                    cds = 0
                    for n in range(i, max(0, i-(3*(win+1))), -1):
                        if record.features[n].type == "CDS" and "pseudo" not in record.features[n].qualifiers:
                            cds += 1
                            if cds <= win and "translation" in record.features[n].qualifiers:
                                orf = SeqRecord(Seq(record.features[n].qualifiers["translation"][0], IUPAC.protein))
                                orf.id = record.features[n].qualifiers["protein_id"][0]
                                SeqIO.write(orf, aa_out, "fasta")
                                upstream = record.features[n].location.start
                                if record.features[n].location.strand == 1:
                                    strand = 'fwd'
                                else:
                                    strand = 'rev'
                                clust[record.features[n].qualifiers["protein_id"][0]] = [strand, int(record.features[n].location.start), int(record.features[n].location.end)]
                    cds = 0
                    if i+1 < len(record.features):
                        for n in range(i+1, min(len(record.features), i+(3*win)), 1):
                            downstream = record.features[n].location.end
                            if record.features[n].type == "CDS" and "pseudo" not in record.features[n].qualifiers:
                                cds += 1
                                #print str(cds)
                                if cds <= win and "translation" in record.features[n].qualifiers:
                                    orf = SeqRecord(Seq(record.features[n].qualifiers["translation"][0], IUPAC.protein))
                                
                                    orf.id = record.features[n].qualifiers["protein_id"][0]
                                    SeqIO.write(orf, aa_out, "fasta")
                                    downstream = record.features[n].location.end
                                    if record.features[n].location.strand == 1:
                                        strand = 'fwd'
                                    else:
                                        strand = 'rev'
                                    clust[record.features[n].qualifiers["protein_id"][0]] = [strand, int(record.features[n].location.start), int(record.features[n].location.end)]
                    else:
                        downstream = record.features[i].location.end
                    nuc_seq = record.seq[upstream:downstream]
                    aa_out.close()
                    print "\t running hmmer"
                    run_hmmer(E, pfam_db)
                    hmmer = open('hmmer.tbl')
                    for hit in hmmer:
                        if hit[0] != "#":
                            clust[hit.split()[0]] = clust[hit.split()[0]]+[hit.split()[3], hit.split()[2], hit.split()[4]]
                    hmmer.close()
                    beginning = []
                    ending = []
                    for k,v in clust.iteritems():
                        beginning.append(v[1])
                        ending.append(v[2])
                    beginning.sort()
                    ending.sort()
                    clust_coords = [beginning[0], ending[-1]]
                    orf_coords = []
                    for num in beginning:
                        for k,v in clust.iteritems():
                            if v[1] == num:
                                protein_acc.writelines(line.qualifiers["protein_id"][0] + ',' + k + ',' + record.id + ',' + ','.join(str(m) for m in v) + '\n')
                                orf_coords.append((v[1], v[2]))
                    
        i += 1
    if "taxonomy" in record.annotations and len(record.annotations["taxonomy"]) > 0:
        phylum = record.annotations["taxonomy"][1]
    else:
        phylum = 'not found'
    if "organism" in record.annotations and len(record.annotations["organism"]) > 0:
        species = record.annotations["organism"].strip()
    else:
        species = "not found"

    nuc = record.seq[clust_coords[0]:clust_coords[1]]
    nuc_acc = record.id
    return p, phylum, species, nuc, clust_coords, orf_coords, nuc_acc

def run_hmmer(E, pfam_db):
    cmd = 'hmmsearch -o hmmer.out --tblout hmmer.tbl -E ' + E + ' --noali ' + pfam_db + ' clust.faa'
    os.system(cmd)
