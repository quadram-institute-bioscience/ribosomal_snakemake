#!/usr/bin/python

def check_for_duplicates(infile):
    import sys
    import re
    from Bio import SeqIO


    handle = open(infile, 'r')
    records = list(SeqIO.parse(handle, "fasta"))

    sequences = []
    outfile = open(sys.argv[1]+"dedupe.fasta", 'w')

    seen_before = []

    def checkers(element, list):
        for l in list:
            if element in l:
                return l


    for rec in records:
        if 'GCF' in rec.id:
            shortid = rec.id.split('.')
            id = shortid[0]
            newid = id
        else:
            newid = rec.id
        if newid in seen_before:
             regex = re.compile(str(newid))
             print("viewed previously", newid, checkers(newid, seen_before))
        else:
             sequences.append(rec)
             seen_before.append(newid)

    SeqIO.write(sequences, outfile, "fasta")

import sys
check_for_duplicates(sys.argv[1])
