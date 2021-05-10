#!/usr/bin/python

from Bio import SeqIO

for rec in records:
...    id = rec.id.split('|')
...    shortid = id[0].split('_')
...    if shortid[1] == ribosomal_proteins[0]:
...         SeqIO.write(rec, out_rplN, "fasta")
...    elif shortid[1] == ribosomal_proteins[1]:
...         SeqIO.write(rec, out_rplP, "fasta")
...    elif shortid[1] == ribosomal_proteins[2]:
...         SeqIO.write(rec, out_rplR, "fasta")
...    elif shortid[1] == ribosomal_proteins[3]:
...         SeqIO.write(rec, out_rplB, "fasta")
...    elif shortid[1] == ribosomal_proteins[4]:
...         SeqIO.write(rec, out_rplV, "fasta")
...    elif shortid[1] == ribosomal_proteins[5]:
...         SeqIO.write(rec, out_rplX, "fasta")
...    elif shortid[1] == ribosomal_proteins[6]:
...         SeqIO.write(rec, out_rplC, 'fasta')
...    elif shortid[1] == ribosomal_proteins[7]:
...         SeqIO.write(rec, out_rplD, 'fasta')
...    elif shortid[1] == ribosomal_proteins[8]:
...         SeqIO.write(rec, out_rplE, 'fasta')
...    elif shortid[1] == ribosomal_proteins[9]:
...         SeqIO.write(rec, out_rplF, 'fasta')
...    elif shortid[1] == ribosomal_proteins[10]:
...         SeqIO.write(rec, out_rpsJ, 'fasta')
...    elif shortid[1] == ribosomal_proteins[11]:
...         SeqIO.write(rec, out_rpsQ, 'fasta')
...    elif shortid[1] == ribosomal_proteins[12]:
...         SeqIO.write(rec, out_rpsS, 'fasta')
...    elif shortid[1] == ribosomal_proteins[13]:
...         SeqIO.write(rec, out_rpsC, 'fasta')
...    elif shortid[1] == ribosomal_proteins[14]:
...         SeqIO.write(rec, out_rpsH, 'fasta')
...    else:
...         print("prob")

