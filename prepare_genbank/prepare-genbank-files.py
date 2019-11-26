from Bio import SeqIO
import os
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

#big_gb = open('ref_seq-hold.gb', 'a')
#big_fasta = open('ref_seq.faa', 'a')
acc_list = open('ref_seq-acc.csv', 'a')
y = 0
z = 0
for i in range(2, 2169):
    print "decompressing files " + str(i)
    os.system('gunzip bacteria.' + str(i) + '.genomic.gbff.gz')
    os.system('gunzip bacteria.' + str(i) + '.1.genomic.fna.gz')
    
    print "indexing files"
    fasta_file = SeqIO.index('bacteria.' + str(i) + '.1.genomic.fna', "fasta")
    gb_file = SeqIO.index('bacteria.' + str(i) + '.genomic.gbff', "genbank")
    
    print "parsing files"
    nuc_hold = open('nuc_hold.gb', 'w')
    for k,record in gb_file.iteritems():
        y += 1
        if record.id in fasta_file:
            if len(fasta_file[record.id].seq) == len(record.seq):
                try:    
                    z += 1
                    record.seq = Seq(str(fasta_file[record.id].seq), IUPAC.ambiguous_dna)
                    SeqIO.write(record, nuc_hold, "genbank")
                    accs = []
                    for feat in record.features:
                        if feat.type == "CDS" and 'protein_id' in feat.qualifiers and 'psuedo' not in feat.qualifiers and 'translation' in feat.qualifiers:
                            #orf = SeqRecord(Seq(feat.qualifiers['translation'][0], IUPAC.protein))
                            #orf.id = feat.qualifiers['protein_id'][0]
                            #SeqIO.write(orf, big_fasta, "fasta")
                            accs.append(feat.qualifiers["protein_id"][0])
                    acc_list.writelines(record.id + ',' + ','.join(accs) + '\n')

                    
                except AttributeError:
                    print "attribute error"
                    continue
        sys.stdout.write('\r%d files processed with %d records and %d over 5,000 bases' % (i, y, z))
        sys.stdout.flush()
    nuc_hold.close()

    print "\nrecompressing file " + str(i)
    #os.system('bgzip nuc_hold.gb')
    os.system('mv nuc_hold.gb bacteria.' + str(i) + '.genomic.gb')
    os.remove('bacteria.' + str(i) + '.genomic.gbff')
    os.remove('bacteria.' + str(i) + '.1.genomic.fna')

gb_file.close()
#fasta_file.close()
acc_list.close()

