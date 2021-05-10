 #!/usr/bin/python


import sys
import re

infile = open(sys.argv[1], 'r')

namesfile = open(sys.argv[2], 'r')

d = {}

for line in namesfile:
    elements = line.rstrip().split('\t')
    taxid = elements[0]
 #    taxid = elements[0][18:28]
  #   print "taxid", taxid
    specname = elements[1]
#    specname = '_'.join(elements[1:5])
  #   print "specname", specname
    try:
        d[taxid].append(specname)
    except:
        d[taxid]=specname

#items = [lin.rstrip() for lin in re.split(r'([^:,\(\) ]+)',line)]
#print(d)

for line in infile:
        items = [lin.rstrip() for lin in re.split(r'([^:,\(\) ]+)',line)]
        for it in items:
           if it in d.keys():
 #                pass
                print(d[it])
           else:
 #               pass
                print(''.join(it))

        

