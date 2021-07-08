#/usr/bin/python

from Bio import SeqIO
import sys

infile = sys.argv[1]

gbk_filename = infile
faa_filename = "{}.faa".format(gbk_filename)
input_handle  = open(gbk_filename, "rt")
output_handle = open(faa_filename, "w")

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
