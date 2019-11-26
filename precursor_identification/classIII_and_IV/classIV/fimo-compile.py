import csv
import math

fh1 = open('fimo.tsv')
next(fh1)

hits = {}

for line in fh1:
    meme = int(line.split()[1].split('-')[1])-1
    if line.split()[2] not in hits:
        scores = [0,0,0,0,0,0,0,0,0]
        scores[meme] = -1*math.log10(float(line.split()[7]))
        hits[line.split()[2]] = scores
    else:
        scores = hits[line.split()[2]]
        if scores[meme] < -1*math.log10(float(line.split()[7])):
            scores[meme] = -1*math.log10(float(line.split()[7]))
            hits[line.strip()[2]] = scores
           
fh1.close()


fh1.close()
fh2 = open('classIV-peps-params.csv')
out = open('classIV-peps-fimo-scores.csv', 'w')
for line in fh2:
    if  "B" not in line.split(',')[7].strip() and "X" not in line.split(',')[7].strip() and "J" not in line.split(',')[7].strip():
        if line.split(',')[0] + '__' + line.split(',')[4] in hits:
            x = [line.split(',')[7].strip()] + hits[line.split(',')[0] + '__' + line.split(',')[4]]
        else:
            x = [line.split(',')[7].strip(),0,0,0,0,0,0,0,0,0]
        csv.writer(out).writerow(x)

fh2.close()
out.close()

