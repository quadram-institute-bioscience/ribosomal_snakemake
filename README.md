# Ribosomal protein tree generation
A workflow to generate ribosomal protein phylogenetic trees, using 15 ribosomal protein sequences.  Either protein or DNA sequences can be used to build the tree, it may be of benefit to use DNA sequences for more closely related genomes, and protein sequences for those that are more divergent.

The original version https://github.com/LCrossman/ribosomal_snakemake by **Lisa Crossman**

This repo has some modifications to make the workflow more portable. Download function is not supported in this version. You are advised to download genome sequences separately using ncbi-genome-download.

# Preparing your data

Suppose you are working on <i>Staphylococcus</i> genus

1.	Collect your ribosomal protein sequences as protein sequence fasta files in 15 separate files into a directory, see `staph_ribosomal_proteins` as an example in this github.
2.	Collect all of your genome sequences as annotated genbank files (either downloaded from genbank NCBI or annotated with Prokka) in another directory, for example `data`.

# Installation and Usage:

On your computer, you are present in a directory with two sub directories: `data` is where you keep your gff/gbk files and `staph_ribosomal_proteins` is for the Staph ribosomal protein.

## Docker
Your computer needs to have [docker](https://docs.docker.com/get-docker/) installed.

```bash
docker run --rm quadram/ribotree ribotree.py --help 
```

If successfully downloaded, you should see the following message:
```
Usage: ribotree.py [OPTIONS] DATA_FOLDER RIBOSOMAL_PROTEIN_FOLDER

Options:
  -c, --config-file TEXT          Config template file  [default:
                                  workflow/config/template_config.yaml]
  -s, --snakefile TEXT            Snakefile  [default: workflow/Snakefile]
  --other-file TEXT               Absolute path to additional ribosomal
                                  protein file
  -t, --threads INTEGER           Number of threads  [default: 8]
  -w, --workflow-dir TEXT         Workflow directory  [default: workflow]
  -o, --outdir TEXT               Output directory  [default: output]
  --tree-builder [iqtree|fasttree]
                                  Tree builder: iqtree or fasttree  [default:
                                  iqtree]
  --protein / --dna               Input type, i.e aminod acid or nucleotide
                                  [default: protein]
  -v, --verbose BOOLEAN           Print verbosity of the execution  [default:
                                  False]
  --dry-run                       Dry run
  --help                          Show this message and exit.
```

```bash
docker run --rm -it -u $(id -u ${USER}):$(id -g ${USER}) -v $PWD:/data quadram/ribotree ribotree.py -t 4  -o test-output ./data ./staph_ribosomal_proteins
```
## Singularity

Your computer needs to have [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html) installed

### Building a Singularity image
```bash
git clone https://github.com/quadram-institute-bioscience/ribosomal_snakemake
sudo singularity build ribotree.sif Singularity
```
### Running the singularity container of the ribotree

Print help
```bash
singularity ribotree.sif ribotree.py --help
```

Execute the pipeline on the 2 folder as an example as mentioned in the Docker section

```
singularity exec ribotree.sif ribotree.py -t 4 -o test-output ./data ./staph_ribosomal_proteins
```
or

```
singularity shell ribotree.sif
ribotree.py -t 4 -o test-output ./data ./staph_ribosomal_proteins 
```
## Conda 

Your computer needs to have [conda](https://conda.io/en/latest/miniconda.html) or [mamba](https://github.com/mamba-org/mamba) installed

```bash
git clone https://github.com/quadram-institute-bioscience/ribosomal_snakemake

# Change to the folder ribosomal_snakemake
cd ribosomal_snakemake
# Install a new conda environment for the pipeline
conda env create -f ribosomal_snakemake/workflow/envs/environment.yaml
# Activate the environment
conda activate rc2_snakemake
# Print help
./ribotree.py --help
# Execute the pipeline on the examples as mentioned in the Docker section
./ribotree.py -t 4 -w workflow -o test-output ./data ./staph_ribosomal_proteins
```
# Pipeline Output

Once the pipeline finishes, you should expect the pipeline output as below

```
test-output/
├── config.yaml
├── logs
│   ├── align.log
│   ├── blast_complete.txt
│   ├── collect_hits.log
│   ├── concatenate.log
│   ├── concatenate_with_previous.log
│   ├── conversion_complete.txt
│   ├── create_tree.log
│   ├── deduplicate.log
│   ├── diamond_run.log
│   ├── extracted_complete.txt
│   └── gunzip_complete.txt
└── results
    ├── 09-07-21.concatenated_ribosomal_proteins_db.fasta
    ├── 09-07-21.concatenated_ribosomal_proteins_db.fasta_2
    ├── 09-07-21.extracted.fasta
    ├── 09-07-21.recovered.fasta
    ├── 09-07-21.updateriboprot.fasta
    ├── 09-07-21.updateriboprot.fastadedupe.fasta
    ├── 09-07-21.updateriboprotdedupe.aln
    ├── 09-07-21.updateriboprotdedupe.aln.treefile
    ├── Allnamesmapoverdatabase.txt
    ├── report.html
    ├── ribosome_db
    └── strains_missing_ribos.txt
```

# TODO
- To make Galaxy wrapper