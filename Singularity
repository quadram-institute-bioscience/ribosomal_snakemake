bootstrap: docker
from: continuumio/miniconda3

%environment
    export PATH=/opt:/opt/conda/envs/rc2_snakemake/bin:$PATH

%files
    workflow /workflow
    ribotree.py /opt/ribotree.py
%post
    export PATH=/opt/conda/bin:$PATH
    apt-get update && apt-get install -y procps
    conda install -c conda-forge mamba
    mamba env create -f /workflow/envs/environment.yaml
