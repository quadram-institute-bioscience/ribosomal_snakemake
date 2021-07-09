#!/usr/bin/env python
from pathlib import Path
import click
from loguru import logger as lg
from snakemake.io import load_configfile
import yaml
from yaml.loader import Loader
import subprocess

def prepare_data_list(folder, ext = ['.gbk','.gb','.genbank','.gff']):
    """
    scan the folder and make a list of file with expected extensions
    """
    genbank_files = Path(folder).glob("*.*")
    results = []
    for genbank_file in genbank_files:
            if genbank_file.suffix in ext:
                lg.info(genbank_file)
                results.append(str(genbank_file))
    return results

def load_config_file(config_file):
    with open(config_file, 'r') as yml_file:
        _config = yaml.load(yml_file, Loader=yaml.FullLoader)
        return(_config)


def run_workflow(snakefile, config_file, threads, dry_run=False):
    cmd = f"snakemake --snakefile {snakefile} --configfile {config_file} --cores {threads}"
    if dry_run:
        cmd += " --dry-run"
    lg.info(f"Started running the workflow: {cmd}")
    
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        lg.critical(e)
        exit(1)

@click.command()
@click.argument('data-folder')
@click.argument('ribosomal-protein-folder')
@click.option('-c','--config-file', default="workflow/config/template_config.yaml", help="Config template file", show_default=True)
@click.option('-s', '--snakefile', default="workflow/Snakefile", help="Snakefile", show_default=True)
@click.option('--other-file', default=None, help="Absolute path to additional ribosomal protein file")
@click.option('-t', '--threads', default=8, help="Number of threads", show_default=True)
@click.option('-w','--workflow-dir', default='workflow', help="Workflow directory", show_default=True)
@click.option('-o', '--outdir', default='output', help="Output directory", show_default=True)
@click.option('--tree-builder',type=click.Choice(['iqtree','fasttree'], case_sensitive=False), default="iqtree", help="Tree builder: iqtree or fasttree", show_default=True)
@click.option('--protein/--dna', default=True, help="Input type, i.e aminod acid or nucleotide", show_default=True)
@click.option('-v','--verbose',default=False, help="Print verbosity of the execution", show_default=True)
@click.option('--dry-run', is_flag=True, default=False, help="Dry run")
def main(data_folder, ribosomal_protein_folder, snakefile, config_file, workflow_dir, outdir, protein, threads, other_file, tree_builder, dry_run, verbose):
    if verbose:
        lg.enable("__main__")
    
    logs_dir = Path(outdir).joinpath("logs")
    lg.info(logs_dir)
    if not logs_dir.exists():
        logs_dir.mkdir(parents=True, exist_ok=True)
    
    input_file = f"{data_folder}/cleannames.txt"
    ribosome_file = f"{data_folder}/atccs.txt"
    this_run_config_file = f"{outdir}/config.yaml"

    data_files = prepare_data_list(data_folder)
    # Prepare a file list of input data
    if len(data_files) > 0:
        with open(input_file, "w") as fh:
            fh.write("\n".join(data_files))

    ribosomal_proteins = prepare_data_list(ribosomal_protein_folder, ext = ['.aa','.fna','.fa'])
    # Prepare a file list of ribosome data
    if len(ribosomal_proteins) > 0:
        with open(ribosome_file, "w") as fh:
            fh.write("\n".join(ribosomal_proteins))

    config = load_config_file(config_file)
    config["cleannames"] = input_file
    config["ribosomefile"] = ribosome_file
    config["outdir"] = outdir
    config["workflow_dir"] = workflow_dir
    config["tree_type"]["options"] = tree_builder
    if not protein:
        config["protein_dna"]["options"] = "dna"
        config["ribo_name_field"]["options"] = 1
    if other_file is not None:
        if Path(other_file).exists():
            config["previous_files"] = other_file
        else:
            lg.enable("__main__")
            lg.critical(f"The additional chromosome protein file {other_file} is not existed.")

    with open(this_run_config_file, "w") as fh:
        yaml.dump(config, fh)

    run_workflow(snakefile=snakefile, config_file = this_run_config_file, threads = threads, dry_run = dry_run)

if __name__ == '__main__':
    lg.disable("__main__")
    main()
