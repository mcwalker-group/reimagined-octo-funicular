fh = open('LanC_like-acc.csv')
classI = open('classI.acc', 'w')
classII = open('classII.acc', 'w')
classIII_IV = open('classIII_IV.acc', 'w')
unknown = open('unknown.acc', 'w')

clusts = {}

for line in fh:
    pfs = []
    stuff = line.strip().split(',')
    for thing in stuff:
        #print thing[0:2]
        if thing[0:2] == "PF":
            pfs.append(thing)
    if stuff[0] in clusts:
        clusts[stuff[0]] = clusts[stuff[0]] + pfs
    else:
        clusts[stuff[0]] = pfs
#print clusts

for k,v in clusts.iteritems():
    if "PF04738.13" in v or "PF14028.6" in v:
        classI.writelines(k + '\n')
    elif "PF13575.6" in v:
        classII.writelines(k + '\n')
    elif "PF05147.13" in v and "PF00069.25" in v:
        classIII_IV.writelines(k + '\n')
    elif "PF05147.13" in v:
        unknown.writelines(k + '\n')

fh.close()
classI.close()
classII.close()
classIII_IV.close()
unknown.close()
