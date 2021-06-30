#!/usr/bin/python
import sys
import collections
import datetime
import os
from Bio import SeqIO
from loguru import logger as lg

def concatenate_diamond_matches(infile, ribo_names_field, outdir):
    handle = open(infile, 'r')
    ribo_names_field = ribo_names_field
    print("ribo_names_field", ribo_names_field)
    d = collections.defaultdict(dict)
    records = list(SeqIO.parse(handle, "fasta"))

    ribosomal_proteins = ['rplN','rplP','rplR','rplB','rplV','rplX','rplC','rplD','rplE','rplF','rpsJ','rpsQ','rpsS','rpsC','rpsH']

##uncomment print statements for debugging, ID changes can cause concatenation to fail

#Loop through annotation and create a dictionary containing each ribosomal sequence per isolate genome

    for rec in records:
        id = rec.id.split('|')
#        print("id to split on |", id)
        shortid = id[0].split('_')
#        print("shortid to split on _", shortid)
        newid = id[1]
#        print("newid", newid)
        rib = shortid[int(ribo_names_field)]
#        print("rib", rib)
        seqy = str(rec.seq)
        if seqy.endswith('*'):
             #remove the stop character sometimes included in protein translation files
             newseq = seqy[:-1]
        else:
             newseq = seqy
        for ribo in ribosomal_proteins:
            if rib == ribo:
                #check for exact match
                try:
                   d[newid][ribo].append(newseq)
                except:
                   d[newid][ribo] = [newseq]

    #create datestamp file
    outname = f"{outdir}/results/{datetime.date.today().strftime('%d-%m-%y')}.concatenated_ribosomal_proteins_db.fasta"
    outfile_strains = open(f"{outdir}/results/strains_missing_ribos.txt", 'a+')
    if os.path.exists(outname):
       #if files have been collected from the annotation and we are now collecting from diamond blast process
       outname = outname+"_2"
    else:
       #mark strains for diamond blast process if they were missing some sequences
       pass

    outfile = open(f"{outname}", 'a+')
    for key, value in d.items():
         #test we have all 15 sequences, otherwise mark genome for diamond blast searches or ignore
         if len(value) == 15:
              joinstring = "{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}".format(''.join(d[key]['rplN']),''.join(d[key]['rplP']),''.join(d[key]['rplR']),''.join(d[key]['rplB']),''.join(d[key]['rplV']),''.join(d[key]['rplX']),''.join(d[key]['rplC']),''.join(d[key]['rplD']),''.join(d[key]['rplE']),''.join(d[key]['rplF']),''.join(d[key]['rpsJ']),''.join(d[key]['rpsQ']),''.join(d[key]['rpsS']),''.join(d[key]['rpsC']),''.join(d[key]['rpsH']))
              outfile.write(">{}\n".format(key))
              outfile.write(joinstring+"\n")
         else:
              print(key, "Missing some ribosomal proteins, will search with Diamond Blast", len(value))
              outfile_strains.write("{}\n".format(key))


if __name__ == "__main__":
     concatenate_diamond_matches(sys.argv[1], sys.argv[2], sys.argv[3])
