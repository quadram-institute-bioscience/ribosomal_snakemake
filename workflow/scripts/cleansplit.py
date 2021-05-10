#!/usr/bin/python

from Bio import SeqIO
import sys

handle = open(sys.argv[1], 'r')
records = list(SeqIO.parse(handle, "fasta"))
ribosomal_proteins = ['rplN','rplP','rplR','rplB','rplV','rplX','rplC','rplD','rplE','rplF','rpsJ','rpsQ','rpsS','rpsC','rpsH']

out_rplN = open("out_rplN.fasta", 'w')
out_rplP = open("out_rplP.fasta", 'w')
out_rplR = open("out_rplR.fasta", 'w')
out_rplB = open("out_rplB.fasta", 'w')
out_rplV = open("out_rplV.fasta", 'w')
out_rplX = open("out_rplX.fasta", 'w')
out_rplC = open("out_rplC.fasta", 'w')
out_rplD = open("out_rplD.fasta", 'w')
out_rplE = open("out_rplE.fasta", 'w')
out_rplF = open("out_rplF.fasta", 'w')
out_rpsJ = open("out_rpsJ.fasta", 'w')
out_rpsQ = open("out_rpsQ.fasta", 'w')
out_rpsS = open("out_rpsS.fasta", 'w')
out_rpsC = open("out_rpsC.fasta", 'w')
out_rpsH = open("out_rpsH.fasta", 'w')






for rec in records:
       id = rec.id.split('|')
       shortid = id[0].split('_')
       if shortid[1] == ribosomal_proteins[0]:
            SeqIO.write(rec, out_rplN, "fasta")
       elif shortid[1] == ribosomal_proteins[1]:
            SeqIO.write(rec, out_rplP, "fasta")
       elif shortid[1] == ribosomal_proteins[2]:
            SeqIO.write(rec, out_rplR, "fasta")
       elif shortid[1] == ribosomal_proteins[3]:
            SeqIO.write(rec, out_rplB, "fasta")
       elif shortid[1] == ribosomal_proteins[4]:
            SeqIO.write(rec, out_rplV, "fasta")
       elif shortid[1] == ribosomal_proteins[5]:
            SeqIO.write(rec, out_rplX, "fasta")
       elif shortid[1] == ribosomal_proteins[6]:
            SeqIO.write(rec, out_rplC, 'fasta')
       elif shortid[1] == ribosomal_proteins[7]:
            SeqIO.write(rec, out_rplD, 'fasta')
       elif shortid[1] == ribosomal_proteins[8]:
            SeqIO.write(rec, out_rplE, 'fasta')
       elif shortid[1] == ribosomal_proteins[9]:
            SeqIO.write(rec, out_rplF, 'fasta')
       elif shortid[1] == ribosomal_proteins[10]:
            SeqIO.write(rec, out_rpsJ, 'fasta')
       elif shortid[1] == ribosomal_proteins[11]:
            SeqIO.write(rec, out_rpsQ, 'fasta')
       elif shortid[1] == ribosomal_proteins[12]:
            SeqIO.write(rec, out_rpsS, 'fasta')
       elif shortid[1] == ribosomal_proteins[13]:
            SeqIO.write(rec, out_rpsC, 'fasta')
       elif shortid[1] == ribosomal_proteins[14]:
            SeqIO.write(rec, out_rpsH, 'fasta')
       else:
            print("prob")

