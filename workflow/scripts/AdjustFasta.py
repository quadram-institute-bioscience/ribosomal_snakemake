#!/usr/bin/python

#Program to create protein fasta files in which the name of the strain is included at the beginning of each protein
#sequence, separated from the locus tag with a | character

import sys

try:
   taxoncode = sys.argv[1]
   inputfield = sys.argv[2]
   idfield = sys.argv[3]
   inp = open(inputfield, 'r')
except:
   print("Usage: python AdjustFasta.py taxoncode inputfield idfield")

outfile = open(inputfield+"_"+taxoncode+".fasta", 'w')

for line in inp:
    if line.startswith('>'):
        details = line.rstrip().split()[1:]
        id = "{}|{}".format(taxoncode, details[int(idfield)])
        outfile.write(">{}\n".format(id))
    else:
        outfile.write(line)
