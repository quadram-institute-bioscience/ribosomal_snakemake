#!/usr/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import sys
import datetime

ribosomal_proteins = ['L14','L16','L18','L2','L22','L24','L3','L4','L5','L6','S10','S17','S19','S3','S8']
ribo_dic = {'L14': 'rplN', 'L16': 'rplP', 'L18': 'rplR', 'L2': 'rplB', 'L22': 'rplV', 'L24': 'rplX', 'L3': 'rplC', 'L4': 'rplD', 'L5': 'rplE', 'L6': 'rplF', 'S10': 'rpsJ', 'S17': 'rpsQ', 'S19': 'rpsS', 'S3': 'rpsC', 'S8': 'rpsH'}

if sys.argv[1].endswith('.gz'):
    filehandle = sys.argv[1][:-3]
    print(filehandle)
else:
    filehandle = sys.argv[1]
handle = open(filehandle, 'r')
params = sys.argv[2]
print("params are", params)
outname = datetime.date.today().strftime("%d-%m-%y")+"extracted.fasta"
records = list(SeqIO.parse(handle, "genbank"))

sequences = []

for rec in records:
    for feature in rec.features:
        if feature.type == "CDS":
            if 'ribosomal protein' in ''.join(feature.qualifiers['product']):
                elements = ''.join(feature.qualifiers['product']).split()
                for rib in ribosomal_proteins:
                    if elements[-1] == rib:
                        print("annotation match")
                        if params == 'dna':
                            newrec = SeqRecord(feature.location.extract(rec.seq), id="{}_{}|{}".format(elements[-1], ribo_dic[elements[-1]], filehandle), description="")
                            sequences.append(newrec)
                        elif params == 'protein':
                            newrec = SeqRecord(feature.location.extract(rec.seq).translate(), id="{}_{}|{}".format(elements[-1], ribo_dic[elements[-1]], filehandle), description="")
                            sequences.append(newrec)
                        else:
                            print("unknown params for DNA or protein")
outfile = open(outname, 'a+')
SeqIO.write(sequences, outfile, "fasta")
