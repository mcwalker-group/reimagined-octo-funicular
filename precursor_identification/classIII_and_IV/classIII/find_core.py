import numpy as np
import re



def leader_motif(fimo_file):
    meme1 = {}
    meme2 = {}
    meme3 = {}
    meme4 = {}
    meme5 = {}
    meme6 = {}
    meme7 = {}
    meme8 = {}
    meme9 = {}

    fimo = open(fimo_file)
    for line in fimo:
        if line[0] != "#":
            if line.split()[1] == "MEME-1" and -1*np.log10(float(line.split()[7])) > 3:
                meme1[line.split()[2]] = (line.split()[3], line.split()[4])
            elif line.split()[1] == "MEME-2" and -1*np.log10(float(line.split()[7])) > 5:
                meme2[line.split()[2]] = (line.split()[3], line.split()[4])
            elif line.split()[1] == "MEME-3" and -1*np.log10(float(line.split()[7])) > 5:
                meme3[line.split()[2]] = (line.split()[3], line.split()[4])
            elif line.split()[1] == "MEME-4" and -1*np.log10(float(line.split()[7])) > 5:
                meme4[line.split()[2]] = (line.split()[3], line.split()[4])
            elif line.split()[1] == "MEME-5" and -1*np.log10(float(line.split()[7])) > 5:
                meme5[line.split()[2]] = (line.split()[3], line.split()[4])
            elif line.split()[1] == "MEME-6" and -1*np.log10(float(line.split()[7])) > 5:
                meme6[line.split()[2]] = (line.split()[3], line.split()[4])
            '''elif line.split()[1] == "MEME-7" and -1*np.log10(float(line.split()[7])) > 5:
                meme7[line.split()[2]] = (line.split()[3], line.split()[4])
            elif line.split()[1] == "MEME-8" and -1*np.log10(float(line.split()[7])) > 5:
                meme8[line.split()[2]] = (line.split()[3], line.split()[4])
            elif line.split()[1] == "MEME-9" and -1*np.log10(float(line.split()[7])) > 5:
                meme9[line.split()[2]] = (line.split()[3], line.split()[4])'''
    fimo.close()
    return [meme1, meme2, meme3]

def cut_pep(pep, lead_motifs):
    acc = pep[0]
    seq = pep[1]

    stops= [0]
    for meme in lead_motifs:
        if acc in meme:
            stops.append(int(meme[acc][0]))
    stops.sort()
    l_stop = stops[-1]

    regex1 = re.search('(.[S|T].{2,7}C)', seq[l_stop:])
    regex2 = re.search('(G[G|A|S])', seq[l_stop:])

    if regex1 and regex2:
        if regex2.end() < regex1.start() and regex2.end()+l_stop > 10:
            leader = seq[0:regex2.end()+l_stop]
            core = seq[regex2.end()+l_stop:]
        elif regex1.start() + l_stop > 10:
            leader = seq[0:regex1.start()+l_stop]
            core = seq[regex1.start()+l_stop:]
        else:
            leader = seq[0:((len(seq)-l_stop)/2)+l_stop]
            core = seq[((len(seq)-l_stop)/2)+l_stop:]
    elif regex1 and not regex2:
        if regex1.start() + l_stop > 10:
            leader = seq[0:regex1.start()+l_stop]
            core = seq[regex1.start()+l_stop:]
        else:
            leader = seq[0:((len(seq)-l_stop)/2)+l_stop]
            core = seq[((len(seq)-l_stop)/2)+l_stop:]
    elif regex2 and not regex1:
        if regex2.end() + l_stop > 10:
            leader = seq[0:regex2.end()+l_stop]
            core = seq[regex2.end()+l_stop:]
        else:
            leader = seq[0:((len(seq)-l_stop)/2)+l_stop]
            core = seq[((len(seq)-l_stop)/2)+l_stop:]            
    else:
        leader = seq[0:((len(seq)-l_stop)/2)+l_stop]
        core = seq[((len(seq)-l_stop)/2)+l_stop:]
    return leader, core


