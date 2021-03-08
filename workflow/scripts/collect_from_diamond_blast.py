#!/usr/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import datetime
import sys
import re

handle = open(sys.argv[1], 'r')
infile_name = sys.argv[1].split('.')
if 'protein' in sys.argv[2]:
   inf = '.'.join(infile_name[2:-1])+".faa"
elif 'dna' in sys.argv[2]:
   inf = '.'.join(infile_name[2:-1])+".ffn"
print("The infile name is", inf)
infile = open(inf, 'r')

keeps = []
d = {}
i = 0
j = []
ribos = ['rplN','rplP','rplR','rplB','rplV','rplX','rplC','rplD','rplE','rplF','rpsJ','rpsQ','rpsS','rpsC','rpsH']

outname = datetime.date.today().strftime("%d-%m-%y")+"recovered.fasta"
print("The outfile name is", outname)

def check_for_more_than_one(handle):
    #checks for more than one diamond blast match
    for line in infile:
        elements = line.rstrip().split()
        j.append(float(elements[-1]))
        maxj = max(j)
        if j.count(maxj) > 1:
            return False
        else:
            return True


d = {}

for ind, line in enumerate(handle):
    elements = line.rstrip().split()
    if any(r for r in ribos if r in elements[1]):
       protid = elements[1]
    else:
       protid = elements[1]+'_'+infile_name[0]
    if check_for_more_than_one == False:
        print("WARNING: there's more than one with the same score!!!")
    else:
        if ind > 0:
            pass
        else:
            try:
               d[elements[0]].append(protid)
            except:
               d[elements[0]] = protid
            keeps.append(elements[0])


#print("this is the keep list", keeps)
            
sequences = []
records = list(SeqIO.parse(infile, "fasta"))
outfile = open(outname, 'a+')
for rec in records:
    try:
        if rec.id in keeps:
            print("Hit identified", rec.id, keeps)
            shortid = rec.id.split('_')
            newid = "{}|{}".format(d[rec.id],shortid[0])
            print("this is the newid", newid)
            new_rec = SeqRecord(rec.seq, id=newid, description="")
            sequences.append(new_rec)
    except:
        pass

print("Number of sequences", len(sequences))
SeqIO.write(sequences, outfile, "fasta")

