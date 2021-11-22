#!/usr/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import datetime
import sys
import re
import os

ext = "ffn"
if sys.argv[2] == "protein":
    ext = "faa"
handle = open(f"{sys.argv[1]}.{ext}", 'r')
infile_name = os.path.basename(sys.argv[1]).split('.')
inf = f"{infile_name[0]}.{ext}"
print("infile name is", inf)
infile = open(f"{sys.argv[1]}.{ext}", 'r')
ribo_names_field = int(sys.argv[3])

keeps = []
d = {}
i = 0
j = []


# outname = f"{datetime.date.today().strftime('%d-%m-%y')}.recovered.fasta"
outname = sys.argv[4]
# print("outname", outname)

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
dicty = {}

for ind, line in enumerate(handle):
    elements = line.rstrip().split()
    # print("elements", elements)
    if check_for_more_than_one == False:
        print("there's more than one with the same score!!!")
    else:
        if ind > 0:
            pass
        else:
            try:
               d[elements[0]].append(elements[1])
            except:
               d[elements[0]] = elements[1]
            keeps.append(elements[0])


print("this is keeps", keeps)

sequences = []
records = list(SeqIO.parse(infile, "fasta"))
outfile = open(outname, 'a+')
for rec in records:
    try:
        if rec.id in keeps:
            print("Hit identified", rec.id, keeps)
            shortid = rec.id.split('_')
            newid = "{}|{}".format(d[rec.id],shortid[int(ribo_names_field)])
            print("this is newid", newid, "here")
            new_rec = SeqRecord(rec.seq, id=newid, description="")
            print("new_rec", new_rec)
            sequences.append(new_rec)
    except:
        pass
print(len(sequences), "len seqs")
SeqIO.write(sequences, outfile, "fasta")
