bootstrap: docker
from: continuumio/miniconda3

%environment
    export PATH=/opt/conda/envs/rc2_snakemake/bin:$PATH

%files
    env_yaml/environment.yaml /environment.yml

%post
    export PATH=/opt/conda/bin:$PATH
    apt-get update && apt-get install -y procps
    conda install -c conda-forge mamba
    mamba env create -f /environment.yml
