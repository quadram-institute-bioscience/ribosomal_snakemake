rule diamond_db:
     input:
          samp=f"{{sample}}",
          miss=f"{config['outdir']}/results/strains_missing_ribos.txt"
     output:
          f"{config['outdir']}/results/ribosome_db/{{sample}}.dmnd"
     shell:
          "diamond makedb --in {input.samp} --db {output}"

rule diamond_run:
     input:
          f"{config['outdir']}/results/strains_missing_ribos.txt"
     output:
          f"{config['outdir']}/logs/blast_complete.txt"
     message: "executing diamond run with {threads} threads on the following files {input}"
     params:
          protein_dna = config['protein_dna']['options']
     log:
          f"{config['outdir']}/logs/diamond_run.log"
     run:
          if os.stat(f"{config['outdir']}/results/strains_missing_ribos.txt").st_size == 0:
               shell(f"touch '{config['outdir']}/logs/blast_complete.txt'")
               shell(f"touch '{config['outdir']}/logs/collect_hits.log'")
               shell(f"touch '{config['outdir']}/logs/concatenate.log'")
          else:
               if params.protein_dna == 'protein':
                    for ribo, name in [(ribo, name) for ribo in SAMPLES for name in [line.rstrip() for line in open('strains_missing_ribos.txt', 'r')]]:
                         shell("diamond blastp --query {0}.faa --db {1} --outfmt 6 --max-target-seqs 1 --out {1}.{0}.out".format(name, ribo))
                    shell(f"touch '{config['outdir']}/logs/blast_complete.txt'")
               else:
                    for ribo, name in [(ribo, name) for ribo in SAMPLES for name in [line.rstrip() for line in open('strains_missing_ribos.txt', 'r')]]:
                         shell("diamond blastx --query {0}.ffn --db {1} --outfmt 6 --max-target-seqs 1 --out {1}.{0}.out".format(name, ribo))
                    shell(f"touch '{config['outdir']}/logs/blast_complete.txt'")

rule collect_hits:
     input:
          infiles=glob.glob("*.out"),
          complete=f"{config['outdir']}/logs/blast_complete.txt"
     output:
          f"{config['outdir']}/results/{DATESTRING['today']}.recovered.fasta"
     params:
          protein_dna = config['protein_dna']['options'],
          ribo_name_field = config['ribo_name_field']['options']
     log:
          log=f"{config['outdir']}/logs/collect_hits.log"
     run:
          print("number of input files", len(input.infiles))
          if len(input.infiles) > 0:
               for inp in input.infiles:
                    print("input is", inp)
                    shell(f"python {config[workflow_dir]}/scripts/collect_from_diamond_blast_nucleotide.py {{inp}} {{params.protein_dna}} {{params.ribo_name_field}}")
          else:
               shell("touch {output}")


rule concatenate:
     input:
          infile=f"{config['outdir']}/results/{DATESTRING['today']}.recovered.fasta"
     output:
          outfile=f"{config['outdir']}/results/{DATESTRING['today']}.concatenated_ribosomal_proteins_db.fasta_2"
     log:
          log=f"{config['outdir']}/logs/concatenate.log"
     params:
          ribo_name_field = config['ribo_name_field']['options']
     shell:
          f"python {config['workflow_dir']}/scripts/ribo_concat_diamond.py {{input.infile}} {{params.ribo_name_field}} {config['outdir']}"
