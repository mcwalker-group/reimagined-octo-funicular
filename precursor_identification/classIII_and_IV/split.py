fh = open('KC_v_L.csv')

fam = {}
next(fh)
for line in fh:
    fam[line.split(',')[0]] = line.strip().split(',')[3]

fh.close()

l = open('classIV-peps.csv', 'w')
kc = open('classIII-peps.csv', 'w')

fh2 = open('classIII_IV-peps.csv')
next(fh2)

for line in fh2:
    if line.split(',')[0] in fam:
        if fam[line.split(',')[0]] == "L":
            l.writelines(line)
        elif fam[line.split(',')[0]] == "KC":
            kc.writelines(line)
fh2.close()
l.close()
kc.close()
