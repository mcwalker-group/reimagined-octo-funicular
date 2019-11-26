import csv

fh1 = open('classI-feature-list.txt')
feats = []
for line in fh1:
    feats.append(line.strip())
print str(len(feats)) + ' features to extract'

fh2 = open('classI-peps-params.csv')
out = open('classI-peps-params-854-feats.csv', 'w')

feats = feats[0:855]

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
