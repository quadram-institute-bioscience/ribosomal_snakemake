configfile: "config.yaml"

import datetime
import glob
import shutil
from pathlib import Path
from os.path import join
datestring = datetime.date.today().strftime("%d-%m-%y")
DATESTRING = {'today': datestring}


include: "download_snake.smk"
include: "blast_missing_seqs.smk"

handle = open("atccs.txt", 'r')
genu = open("genus.txt", 'r')
GENUS = [liner.rstrip() for liner in genu]
SAMPLES = [lin.rstrip() for lin in handle]

print(DATESTRING, SAMPLES, GENUS)

rule all:
    input:
        "logs/cleannamecheck.txt",
        "logs/completed.txt",
        "logs/newroot.txt",
        "logs/gunzip_complete.txt",
        "logs/conversion_complete.txt",
        "Allnamesmapoverdatabase.txt",
        atccs=expand("{sample}.dmnd", sample=SAMPLES),
        missing="strains_missing_ribos.txt",
        conc=DATESTRING['today']+"concatenated_ribosomal_proteins_db.fasta",
        extr="logs/extracted_complete.txt",
        subtest=DATESTRING['today']+"concatenated_ribosomal_proteins_db.fasta_2",
        dmnd="logs/blast_complete.txt",
        report="report.html"


def glob_files():
      suffixes = ['*.gbff', '*gbf', '*.gbk.gbk', '*.gbk', '*.gb', '*.genbank']
      files_glob = []
      for suf in suffixes:
          files_glob.append([Path(fa).name for fa in glob.glob(suf)])
      if 'no' in config['download_genbank']['options']:
          shell("touch 'logs/gunzip_complete.txt'")
      else:
          pass
      return files_glob

rule check_for_cleannames:
      input:
           files = glob.glob('*.txt')
      output:
           "logs/cleannamecheck.txt"
      params:
           required = config['download_genbank']['options']
      run:
           shell("touch 'cleannames.txt'")
           if 'yes' in params.required:
               shell("touch 'cleannames.txt'")
           else:
               if 'cleannames.txt' in input.files:
                   pass
               else:
                   print("Please create a file called cleannames.txt with the genbank files you want to include in the analysis one filename per line")
           shell("touch 'logs/cleannamecheck.txt'")

rule convert_nucl_protein:
      input:
           files = glob_files(),
           gunz = "logs/gunzip_complete.txt",
           clean = "logs/cleannamecheck.txt"
      output:
           "logs/conversion_complete.txt"
      params:
           type=config['protein_dna']['options'],
           required=config['download_genbank']['options']
      threads:
           config['threads']
      run:
           if 'yes' in params.required:
               outfile2 = open("cleannames.txt", 'w')
           for file in input.files:
               print("this is file", file)
               if 'yes' in params.required:
                   outfile2.write(Path(file).name+"\n")
               if 'protein' in params.type:
                   shell("python gbk2faa.py {file}")
               elif 'dna' in params.type:
                   shell("python gbk2ffn.py {file}")
           shell("touch 'logs/conversion_complete.txt'")

rule create_mapover:
       input:
            "logs/conversion_complete.txt"
       output:
            "Allnamesmapoverdatabase.txt"
       params:
            required=config['download_genbank']['options']
       run:
            if 'yes' in params.required:	  
                shell("for file in GCF*gbff; do echo $file; grep DEFINITION $file; done > whichsequenceiswhich.txt")
                shell("python makeconsdatabasemapping.py > Allnamesmapoverdatabase.txt")
            else:
                shell("touch 'Allnamesmapoverdatabase.txt'")


rule extract_from_gbk:
    input:
        check=glob.glob("cleannames.txt"),
        mapov="Allnamesmapoverdatabase.txt"
    output:
        DATESTRING['today']+"extracted.fasta",
        "logs/extracted_complete.txt"
    params:
        config['protein_dna']['options']
    run:
        inputs = [lin.rstrip() for lin in open(''.join(input.check),'r')]
        for cleanname in inputs:
            shell("python extract_ribo_seqs_from_gbk.py {cleanname} {params}")
        shell("touch 'logs/extracted_complete.txt'")


rule check_for_missing_seqs:
    input:
        extracted=DATESTRING['today']+"extracted.fasta",
        comple="logs/extracted_complete.txt"
    output:
        DATESTRING['today']+"concatenated_ribosomal_proteins_db.fasta",
        "strains_missing_ribos.txt"
    params:
        protein_dna = config['protein_dna']['options'],
        ribo_name_field = config['ribo_name_field']['options']
    run:
        shell(f"python ribo_concat_diamond.py {input.extracted} {params.ribo_name_field}")


rule concatenate_with_previous:
    input:
         cats = config["previous_files"],
         newput = DATESTRING['today']+"concatenated_ribosomal_proteins_db.fasta",
         nextput = DATESTRING['today']+"concatenated_ribosomal_proteins_db.fasta_2"
    output:
         DATESTRING['today']+".updateriboprot.fasta"
    log:
         log="logs/concatenate_with_previous.log"
    shell:
         "(cat {input.cats} {input.newput} {input.nextput} > {output}) 2>>log"


rule deduplicate:
     input:
         DATESTRING['today']+".updateriboprot.fasta"
     output:
         dedupe = DATESTRING['today']+".updateriboprot.fastadedupe.fasta",
     log:
         log="logs/deduplicate.log"
     shell:
         "python check_not_duplicatedIDs.py {input}"

rule align:
     input:
         DATESTRING['today']+".updateriboprot.fastadedupe.fasta"
     output:
         DATESTRING['today']+".updateriboprotdedupe.aln"
     threads:
         config['threads']
     log:
         log="logs/align.log"
     run:
         shell("(mafft --retree 3 --maxiterate 3 --thread {threads} {input} > {output}) 2>>log")

rule create_tree:
      input:
          DATESTRING['today']+".updateriboprotdedupe.aln"
      output:
          DATESTRING['today']+".updateriboprotdedupe.aln.treefile"
      params:
          tree_type = config['tree_type']['options']
      threads:
          config['threads']
      log:
          "logs/create_tree.log"
      run:
          if params.tree_type == 'iqtree':
               shell("(iqtree -s {input} -m LG -bb 1000 -alrt 1000 -nt {threads} > {output}) 2>>log")
          else:
               shell("(fasttree < {input} > {output}) 2>>log")

rule report:
      input:
          DATESTRING['today']+".updateriboprotdedupe.aln.treefile"
      output:
          "report.html"
      run:
         from snakemake.utils import report
         with open(input[0]) as aln:
            n_strains = sum(1 for a in aln if a.startswith('>'))
         report("""A workflow initially formulated for downloading Non-aureus Staph genomes, processing to extract ribosomal protein sequences, deduplicating, aligning and creating a phylogenetic tree""",output[0], T1=input[0])
