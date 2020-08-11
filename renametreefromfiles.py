 #!/usr/bin/python


import sys
import re

infile = open(sys.argv[1], 'r')

namesfile = open(sys.argv[2], 'r')

d = {}

for line in namesfile:
    elements = line.rstrip().split('\t')
    taxid = elements[1]
 #    taxid = elements[0][18:28]
    specname = elements[0]
#    specname = '_'.join(elements[1:5])
    try:
        d[taxid].append(specname)
    except:
        d[taxid]=specname

for line in infile:
        items = [lin.rstrip() for lin in re.split(r'([^:,\(\) ]+)',line)]
        for it in items:
           if it in d.keys():
                print(d[it])
           else:
                print(''.join(it))
