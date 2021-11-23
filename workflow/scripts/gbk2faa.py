#/usr/bin/python

from Bio import SeqIO
import sys
import os

outdir = sys.argv[2]

handle_name = os.path.basename(sys.argv[1]).split('.')[0]
outfile = os.path.join(outdir,f"{handle_name}.faa")

input_handle = open(sys.argv[1], 'rt')
output_handle = open(outfile, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record {}".format(seq_record.id))
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            try:
               assert len(seq_feature.qualifiers['translation'])==1
               output_handle.write(">{} from {}\n{}\n".format(
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_record.name,
                   seq_feature.qualifiers['translation'][0]))
            except:
                print("missing value", seq_feature)

output_handle.close()
input_handle.close()
