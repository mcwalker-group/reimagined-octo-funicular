fh = open('lanC-like.acc')

lanc = []

for line in fh:
    lanc.append(line.strip())
fh.close()

x = len(lanc)
fh2 = open('/mnt/d/bacteria_refseq/ref_seq-acc.csv')
out = open('lanC-like_gb_ids.csv', 'w')
y = 0
for line in fh2:
    y += 1
    if len(set(lanc) & set(line.strip().split(','))) > 0:
        out.writelines(line)
        x = x-1
    print str(x) + ' left to find after searching ' + str(y) + ' records'
