#!/usr/bin/python

def concatenate_diamond_matches(infile):
    import collections
    import datetime
    import os
    from Bio import SeqIO

    handle = open(infile, 'r')
    d = collections.defaultdict(dict)
    records = list(SeqIO.parse(handle, "fasta"))

    ribosomal_proteins = ['rplN','rplP','rplR','rplB','rplV','rplX','rplC','rplD','rplE','rplF','rpsJ','rpsQ','rpsS','rpsC','rpsH']

    for rec in records:
        id = rec.id.split('|')
#        print("id to split on |", id)
        shortid = id[0].split('_')
#        print("shortid to split on _", shortid)
        newid = id[1]
#        print("newid", newid)
        rib = shortid[1]
#        print("rib", rib)
        seqy = str(rec.seq)
        if seqy.endswith('*'):
             newseq = seqy[:-1]
        else:
             newseq = seqy
        for ribo in ribosomal_proteins:
            if rib == ribo:
                try:
                   d[newid][ribo].append(newseq)
                except:
                   d[newid][ribo] = [newseq]

    outname = datetime.date.today().strftime("%d-%m-%y")+"concatenated_ribosomal_proteins_db.fasta"
    if os.path.exists(outname):
       outname = outname+"_2"
    else:
       outfile_strains = open("strains_missing_ribos.txt", 'a+')
    outfile = open(outname, 'a+')
    for key, value in d.items():
         if len(value) == 15:
              joinstring = "{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}".format(''.join(d[key]['rplN']),''.join(d[key]['rplP']),''.join(d[key]['rplR']),''.join(d[key]['rplB']),''.join(d[key]['rplV']),''.join(d[key]['rplX']),''.join(d[key]['rplC']),''.join(d[key]['rplD']),''.join(d[key]['rplE']),''.join(d[key]['rplF']),''.join(d[key]['rpsJ']),''.join(d[key]['rpsQ']),''.join(d[key]['rpsS']),''.join(d[key]['rpsC']),''.join(d[key]['rpsH']))
              outfile.write(">{}\n".format(key))
              outfile.write(joinstring+"\n")
         else:
              print(key, "Missing some ribosomal proteins, will search with Diamond Blast", len(value))
              outfile_strains.write("{}\n".format(key))

import sys
concatenate_diamond_matches(sys.argv[1])
