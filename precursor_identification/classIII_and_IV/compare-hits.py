fh1 = open('III_IV.acc')

acc = []

for line in fh1:
    acc.append(line.strip())
fh1.close()

fh2 = open('lanKC.tbl')
kc = {}
for line in fh2:
    if line[0] != "#":
        kc[line.split()[0]] = line.split()[4]
fh2.close()

fh3 = open('lanL.tbl')

l = {}

for line in fh3:
    if line[0] != "#":
        l[line.split()[0]] = line.split()[4]
fh3.close()

for thing in acc:
    if thing in kc and thing in l:
        print thing + ',' + kc[thing] + ',' + l[thing]
    elif thing in kc and thing not in l:
        print thing + ',' + kc[thing] + ',100'
    elif thing not in kc and thing in l:
        print thing + ',100' + ',' + l[thing]
    else:
        print thing + ',100,100'



