#!/usr/bin/python3


handle = open("whichsequenceiswhich.txt", 'r')
file = []
species = []
i=0


for line in handle:
    if line.startswith('GCF'):
        filename = line.rstrip()
        fil = filename.split('_')
        shortname = "GCF_"+fil[1]
        file.append(filename)
        i=0
    elif line.startswith('DEFINITION'):
        i+=1
        if i > 1:
            pass
        else:
            elements = line.split()
            species.append("{}_{}".format(shortname,elements[2]))
    else:
        print("odd error")

    
mapover = zip(file, species)

for m in mapover:
    print ("{}".format('\t'.join(m)))
