#code for peptide processing

import re

def get_peps(begin, nuc, direction, pep_max, pep_min):

    starts = []
    peps = []
    for i in range(0,len(nuc),3):
        codon = nuc[i:i+3]
        if codon == "ATG" or codon == "GTG" or codon == "TTG":
            starts.append(i)

    for x in range(0, len(starts)):
        regex = re.search('^((ATG|GTG|TTG)(\w\w\w){0,%d}?(TAA|TGA|TAG))' % pep_max, str(nuc[starts[x]:]))
        if regex:
            if len(nuc[regex.start()+starts[x]:regex.end()+starts[x]].translate()) >= pep_min and len(nuc[regex.start()+starts[x]:regex.end()+starts[x]].translate()) <= pep_max:
                if direction == 'f':
                    peps.append((begin+starts[x], begin+starts[x]+regex.end(), str(nuc[regex.start()+starts[x]:regex.end()+starts[x]-3].translate())))
                if direction == 'r':
                     peps.append((begin-starts[x]-regex.end(), begin-starts[x], str(nuc[regex.start()+starts[x]:regex.end()+starts[x]-3].translate())))
                
    return peps

def remove_redundancy(peps, direction):
    longest = {}
    if direction == 'f':
        for pre in peps:
            if pre[1] not in longest:
                longest[pre[1]] = pre
            elif len(longest[pre[1]][2]) < len(pre[2]):
                longest[pre[1]] = pre
    if direction == 'r':
        for pre in peps:
            if pre[0] not in longest:
                longest[pre[0]] = pre
            elif len(longest[pre[0]][2]) < len(pre[2]):
                longest[pre[0]] = pre
    return longest.values()

def intergenic(f, r, coords):
    inter_f = []
    inter_r = []
    for pep in f:
        trash = 'no'
        for orf in coords:
            if pep[0] >= orf[0] and pep[1] <= orf[1]:  #peptide ORF is inside other orf
                trash = 'yes'
            if pep[0] == orf[0] and pep[1] == orf[1]:   #peptide orf is other orf
                trash = 'no'
            #elif pep[0] < orf[0] and pep[1] == orf[1]:     #peptide orf starts before other orf but ends in the same spot
            #    trash = 'no'
            #elif pep[0] < orf[0] and pep[1] > orf[0] and pep[1] - 12 > orf[0]:   #peptide orf starts 5' to other org and overlaps by more than 12 bases
            #    trash = 'yes'
            #elif pep[0] < orf[1] and pep[1] > orf[1] and pep[0] + 12 < orf[1]:    #peptide orf ends 3' to other org but starts more than 12 bases inside other orf
            #    trash = 'yes'
        if trash == 'no':
            inter_f.append(pep)
    for pep in r:
        trash = 'no'
        for orf in coords:
            if pep[0] >= orf[0] and pep[1] <= orf[1]:
                trash = 'yes'
            if pep[0] == orf[0] and pep[1] == orf[1]:   #peptide orf is other orf
                trash = 'no'
            #elif pep[1] > orf[1] and pep[0] == orf[0]:     #peptide orf starts before other orf but ends in the same spot
            #    trash = 'no'
            #elif pep[0] < orf[0] and pep[1] > orf[0] and pep[1] - 12 > orf[0]:
            #    trash = 'yes'
            #elif pep[0] < orf[1] and pep[1] > orf[1] and pep[0] + 12 < orf[1]:
            #    trash = 'yes'
        if trash == 'no':
            inter_r.append(pep)    
    return inter_f, inter_r
