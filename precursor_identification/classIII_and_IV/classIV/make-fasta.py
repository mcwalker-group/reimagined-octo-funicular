fh = open('classIV-precursors-dereplicated.csv')
fh.next()
for line in fh:
    print ">"+line.split(',')[0] + '__' + line.split(',')[4] + '\n' + line.split(',')[9].strip()
