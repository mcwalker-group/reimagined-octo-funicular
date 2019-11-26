import csv

fh1 = open('classIV-feat-list.txt')
feats = []
for line in fh1:
    feats.append(line.strip())
print str(len(feats)) + ' features to extract'

fh2 = open('classIV-peps-params.csv')
out = open('classIV-peps-params-1007-feats.csv', 'w')

feats = feats[0:1008]

locs = []
for line in fh2:
    stuff = line.strip().split(',')
    if stuff[0] == 'query_acc':
        for thing in feats:
            locs.append(stuff.index(thing))
        new = stuff[0:2] + feats
        csv.writer(out).writerow(new)
    else:
        new = stuff[0:2]
        for num in locs:
            new.append(stuff[num])
        csv.writer(out).writerow(new)
fh2.close()
out.close()
