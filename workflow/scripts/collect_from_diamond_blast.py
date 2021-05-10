#!/usr/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter
from toolz import compose
import csv
import datetime
import sys
import re

handle = open(sys.argv[1], 'r')
reader = csv.reader(open(sys.argv[1]), delimiter='\t')
infile_name = sys.argv[1].split('.')
if 'protein' in sys.argv[2]:
   inf = '.'.join(infile_name[2:-1])+".faa"
elif 'dna' in sys.argv[2]:
   inf = '.'.join(infile_name[2:-1])+".ffn"
print("infile name is", inf)
infile = open(inf, 'r')

keeps = []
d = {}
i = 0
j = []

outname = datetime.date.today().strftime("%d-%m-%y")+"recovered.fasta"
print("outname", outname)

def check_for_more_than_one(handle):
    for line in infile:
        elements = line.rstrip().split()
        j.append(float(elements[-1]))
        maxj = max(j)
        if j.count(maxj) > 1:
            return False
        else:
            return True


d = {}
#reader = csv.reader(open(handle), delimiter='\t')

for i, line in enumerate(sorted(reader, key=compose(float, itemgetter(2)), reverse=True)):
    print(line)
   # elements = line.rstrip().split()
    if check_for_more_than_one == False:
        print("there's more than one with the same score!!!")
    else:
        if i > 0:
            pass
        else:
            try:
               d[line[0]].append(line[1])
            except:
               d[line[0]] = line[1]
            keeps.append(line[0])


print("this is keeps", keeps)
            
sequences = []
records = list(SeqIO.parse(infile, "fasta"))
outfile = open(outname, 'a+')
for rec in records:
    try:
        if rec.id in keeps:
            print("Hit identified", rec.id, keeps)
            shortid = rec.id.split('_')
            newid = "{}|{}".format(d[rec.id],shortid[0])
            print("this is newid", newid, "here")
            new_rec = SeqRecord(rec.seq, id=newid, description="")
            print("new_rec", new_rec)
            sequences.append(new_rec)
    except:
        pass

print("len seqeunces", len(sequences))
SeqIO.write(sequences, outfile, "fasta")

