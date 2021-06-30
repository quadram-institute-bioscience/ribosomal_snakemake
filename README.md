# ribosomal_snakemake

# Ribosomal protein tree generation

A workflow to generate ribosomal protein phylogenetic trees, using 15 ribosomal protein sequences.  Either protein or DNA sequences can be used to build the tree, it may be of benefit to use DNA sequences for more closely related genomes, and protein sequences for those that are more divergent.

<img src="https://github.com/LCrossman/ribosomal_snakemake/blob/main/RC_snakemake2.jpg" width="450" height="550">

# Requirements:
snakemake>=3.5<br>
python>=3.3<br>
biopython>=1.77<br>
toolz>=0.11.1<br>
diamond BLAST version>=0.9.32<br>
mafft>=7.429<br>
fasttree>=2.1.10 OR:<br>
iqtree>=1.6.11<br>
ncbi-genome-download>=0.3.0

Python scripts from this github folder

An environment.yaml is included for conda in env_yaml directory if required.

# Installation <i>via</i> conda:
Use conda to install snakemake to a linux or Mac OSX environment:
Install conda:<br>
Install conda here https://conda.io/en/latest/miniconda.html<br>
Install dependencies <i>via</i> conda 

Install snakemake:<br>

conda install -c conda-forge mamba<br>
mamba create -c conda-forge -c bioconda -n snakemake snakemake<br>
<br>
conda activate snakemake<br>
snakemake --help<br>

Install the ribosomal protein tree workflow from github:<br>
git clone XXX<br>
cd XXX<br>
cd testfiles<br>
snakemake -n<br>

# Usage
Configure workflow:<br>
Configure the workflow according to your needs by editing the file config.yaml. For use without any sequence database downloads, you need to provide the files and pathnames in cleannames.txt and atccs.txt.

1.	Gather your ribosomal protein sequences as protein sequence fasta files in 15 separate files.  These will be used for either protein or DNA sequence trees as per your choice laid out in the config.yaml.  They should be named within the file as L14_rplN, L16_rplP, L18_rplR, L2_rplB, L22_rplV, L24_rplX, L3_rplC, L4_rplD, L5_rplE, L6_rplF, S10_rpsJ, S17_rpsQ, S19_rpsS, S3_rpsC and S8_rpsH, respectively and 1 should be placed in the ribo_names_field in the config.yaml.  If they are named with the rpl or rpS letter first (such as rplN_L14) 0 should be placed in the config.yaml file under ribo_name_field.  <i>Staphylococcus</i> sequences are included in the github folder and for examples.
2.	Gather all of your genome sequences as annotated genbank files in the same directory. 
3.	Create a file, “genus.txt” containing the name of the genus and species you want to build a tree for, one per line.  Note that the genus or species name should be enclosed with quotes.
4.	 If you do not want to download any genomes from the sequence databases, create a file, “cleannames.txt” containing the name of each genome file you want to include in your tree, one per line, and place “no” in the config.yaml under download_genbank options.  To use a mixture of download and provided genomes place “yes” in the config.yaml and do not provide a cleannames.txt file.
5.	To update any previous trees, collect any files generated as concatenated deduplicated fasta files that you want to update, and edit appropriately the config.yaml file with the filenames.  
Add Genbank files, ribosomal protein sequence files and a list of the files:
If you do not want to download any genomes, add the names of all provided genbank format files to cleannames.txt on a one-line per file basis.  The suffix of the genbank files should end in .gbff (following the genbank download nomenclature), .gbk, .gb or .genbank. Add the names of your ribosomal sequence files on a one-line per file basis to atccs.txt.  An easy way to generate these files is with:
ls *gbff > cleannames.txt
<br>
Example files are contained in example_files folder with files for a testrun in the testfiles folder
# Execute workflow:
Test your configuration by performing a dry-run via
snakemake -n<br>
Execute the workflow locally via<br>
snakemake --cores $N<br>
using $N cores<br>

or run it in a cluster environment such as<br>
snakemake --use-conda --cluster qsub --jobs 100<br>

Further information on running snakemake in a cluster environment can be found on the snakemake website<br>
