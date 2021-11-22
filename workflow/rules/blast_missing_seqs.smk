import os

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
        strains_missing_ribos = f"{config['outdir']}/results/strains_missing_ribos.txt"
        # converted = f"{config['outdir']}/logs/conversion_complete.txt"
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
            shell(f"mkdir -p {config['outdir']}/diamond")
            _strains_missing_ribos = [line.rstrip() for line in open(input.strains_missing_ribos, 'r')]
            if params.protein_dna == 'protein':
                for ribo, name in [(ribo, name) for ribo in SAMPLES for name in _strains_missing_ribos]:
                    ribo_name = os.path.basename(ribo).split(".")[0]
                    shell(f"diamond blastp --query {config['converted_nuc']}/{name}.faa --db {ribo} --outfmt 6 --max-target-seqs 1 --out {config['outdir']}/diamond/{name}.{ribo_name}.out")
                shell(f"touch '{config['outdir']}/logs/blast_complete.txt'")
            else:
                for ribo, name in [(ribo, name) for ribo in SAMPLES for name in _strains_missing_ribos]:
                    ribo_name = os.path.basename(ribo).split(".")[0]
                    shell(f"diamond blastx --query {config['converted_nuc']}/{name}.ffn --db {ribo} --outfmt 6 --max-target-seqs 1 --out {config['outdir']}/diamond/{name}.{ribo_name}.out")
                shell(f"touch '{config['outdir']}/logs/blast_complete.txt'")

rule collect_hits:
    input:
        infiles=glob.glob(f"{config['outdir']}/diamond/*.out"),
        complete=f"{config['outdir']}/logs/blast_complete.txt"
    output:
        file_name = f"{config['outdir']}/results/{DATESTRING['today']}.recovered.fasta"
    params:
        protein_dna = config['protein_dna']['options'],
        ribo_name_field = config['ribo_name_field']['options']
    log:
        log=f"{config['outdir']}/logs/collect_hits.log"
    run:
        # print("number of input files", len(input.infiles))
        if len(input.infiles) > 0:
            for _input in input.infiles:
                _input_file = os.path.basename(_input).split(".")[0]
                shell(f"python {config['workflow_dir']}/scripts/collect_from_diamond_blast_nucleotide.py {config['converted_nuc']}/{_input_file} {{params.protein_dna}} {{params.ribo_name_field}} {{output.file_name}}")
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
