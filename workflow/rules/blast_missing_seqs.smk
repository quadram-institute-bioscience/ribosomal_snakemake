rule diamond_db:
    input:
         samp=f"{{sample}}",
	 miss="strains_missing_ribos.txt"
    output:
         f"{{sample}}.dmnd"
    shell:
         "diamond makedb --in {input.samp} --db {output}"

rule diamond_run:
     input:
         "strains_missing_ribos.txt"
     output:
         "logs/blast_complete.txt"
     message: "executing diamond run with {threads} threads on the following files {input}"
     params:
          protein_dna = config['protein_dna']['options']
     log:
          "logs/diamond_run.log"
     run:
          if os.stat('strains_missing_ribos.txt').st_size == 0:
               shell("touch 'logs/blast_complete.txt'")
               shell("touch 'logs/collect_hits.log'")
               shell("touch 'logs/concatenate.log'")
          else:
              if params.protein_dna == 'protein':
                   for ribo, name in [(ribo, name) for ribo in SAMPLES for name in [line.rstrip() for line in open('strains_missing_ribos.txt', 'r')]]:
                       shell("diamond blastp --query {0}.faa --db {1} --outfmt 6 --max-target-seqs 1 --out {1}.{0}.out".format(name, ribo))
                   shell("touch 'logs/blast_complete.txt'")
              else:
                   for ribo, name in [(ribo, name) for ribo in SAMPLES for name in [line.rstrip() for line in open('strains_missing_ribos.txt', 'r')]]:
                       shell("diamond blastx --query {0}.ffn --db {1} --outfmt 6 --max-target-seqs 1 --out {1}.{0}.out".format(name, ribo))
                   shell("touch 'logs/blast_complete.txt'")

rule collect_hits:
    input:
         infiles=glob.glob("*.out"),
	 complete="logs/blast_complete.txt"
    output:
         DATESTRING['today']+"recovered.fasta"
    params:
         protein_dna = config['protein_dna']['options']
    log:
         log="logs/collect_hits.log"
    run:
         print("len", len(input.infiles))
         if len(input.infiles) > 0:
             for inp in input.infiles:
                 print("inp is", inp)
                 if params.protein_dna == 'protein':
                     shell(f"python collect_from_diamond_blast.py {inp} {params.protein_dna}")
                 else:
                     shell(f"python collect_from_diamond_blast_nucleotide.py {inp} {params.protein_dna}")
         else:
             shell("touch {output}")
		   

rule concatenate:
    input:
         infile=DATESTRING['today']+"recovered.fasta"
    output:
         outfile=DATESTRING['today']+"concatenated_ribosomal_proteins_db.fasta_2"
    log:
         log="logs/concatenate.log"
    shell:
         "python ribo_concat_diamond.py {input.infile}"
