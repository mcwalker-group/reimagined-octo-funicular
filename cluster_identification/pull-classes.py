fh1 = open('classII.acc')
acc = []
for line in fh1:
    acc.append(line.strip())
fh1.close()

fh2 = open('lanC_like-acc.csv')
for line in fh2:
    if line.split(',')[0] in acc:
        print line.strip()