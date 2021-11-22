#!/usr/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import sys
import datetime
from pathlib import Path
from loguru import logger as lg

ribosomal_proteins = ['L14','L16','L18','L2','L22','L24','L3','L4','L5','L6','S10','S17','S19','S3','S8']
ribo_dic = {'L14': 'rplN', 'L16': 'rplP', 'L18': 'rplR', 'L2': 'rplB', 'L22': 'rplV', 'L24': 'rplX', 'L3': 'rplC', 'L4': 'rplD', 'L5': 'rplE', 'L6': 'rplF', 'S10': 'rpsJ', 'S17': 'rpsQ', 'S19': 'rpsS', 'S3': 'rpsC', 'S8': 'rpsH'}

if sys.argv[1].endswith('.gz'):
    filehandle = sys.argv[1][:-3]
    print(filehandle)
else:
    filehandle = sys.argv[1]

handle = open(filehandle, 'r')

params = sys.argv[2]
outdir = sys.argv[3]
# lg.info(f"parameters for protein or dna are: {params}")
outname = f"{datetime.date.today().strftime('%d-%m-%y')}.extracted.fasta"
records = list(SeqIO.parse(handle, "genbank"))

sequences = []

#Loop through annotation to look for an annotation match
#Save correct matches to file
filename = Path(filehandle).stem
for rec in records:
    for feature in rec.features:
        if feature.type == "CDS":
            if 'ribosomal protein' in ''.join(feature.qualifiers['product']):
                elements = ''.join(feature.qualifiers['product']).split()
                for rib in ribosomal_proteins:
                    if elements[-1] == rib:
                        # lg.info(f"annotation match: {rib[-1]}")
                        if params.lower() == 'dna':
                            #DNA parameters so extract the DNA sequence and save to a file with the appropriate ID
                            newrec = SeqRecord(feature.location.extract(rec.seq), id="{}_{}|{}".format(
                                elements[-1], ribo_dic[elements[-1]], filename), description="")
                            sequences.append(newrec)
                        elif params.lower() == 'protein':
                            #protein parameters so extract the protein sequence, (by translating the DNA sequence) and save to file with appropriate ID
                            newrec = SeqRecord(feature.location.extract(rec.seq).translate(), id="{}_{}|{}".format(
                                ribo_dic[elements[-1]], elements[-1], filename), description="")
                            sequences.append(newrec)
                        else:
                            lg.info(f"unknown parameters for DNA or protein, please choose either dna or protein")

lg.info(f"{filehandle} has {len(sequences)} squences extracted!")

outfile = open(f"{outdir}/results/{outname}", 'a+')
#save chosen sequences to file with date stamp

SeqIO.write(sequences, outfile, "fasta")
