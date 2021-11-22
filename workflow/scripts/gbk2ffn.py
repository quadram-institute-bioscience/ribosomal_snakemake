#!/usr/bin/python


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import sys
import os


outdir = sys.argv[2]
handle = open(sys.argv[1], 'rt')
handle_name = os.path.basename(sys.argv[1]).split('.')[0]
outfile = open(os.path.join(outdir,f"{handle_name}.ffn"), 'w')

sequences = []
records = list(SeqIO.parse(handle, "genbank"))

for rec in records:
    # print("Dealing with genbank record {}".format(rec.id))
    for feature in rec.features:
        if feature.type == "CDS":
            if feature.location.strand == -1:
                ffn = rec.seq[feature.location.start:feature.location.end]
                ffn.r = ffn[::-1]
                newid = ''.join(feature.qualifiers['locus_tag'])
                newdescription = ''.join(feature.qualifiers['product'])
                sequences.append(SeqRecord(ffn.r,id=newid,description=newdescription))
            elif feature.location.strand == 1:
                ffn = rec.seq[feature.location.start:feature.location.end]
                newid = ''.join(feature.qualifiers['locus_tag'])
                newdescription = ''.join(feature.qualifiers['product'])
                sequences.append(SeqRecord(ffn,id=newid,description=newdescription))

SeqIO.write(sequences, outfile, "fasta")
