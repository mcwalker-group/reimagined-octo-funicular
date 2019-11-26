import csv
import math

fh1 = open('leader-hmm.tbl')
acc = []
for line in fh1:
    if line[0] != "#":
        if line.split()[2] not in acc:
            acc.append(line.split()[2])

fh1.close()

fh2 = open('classII-peps-params.csv')
out = open('classII-peps-HMM.csv', 'w')
i = 0
for line in fh2:
    i += 1
    if line.split(',')[0] + "__" + line.split(",")[4] in acc:
        new = [i, 1]
        csv.writer(out).writerow(new)
    else:
        new = [i, 0]
        csv.writer(out).writerow(new)
fh2.close()
out.close()

