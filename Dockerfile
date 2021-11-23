FROM ubuntu:18.04

LABEL org.opencontainers.image.authors="thanh.le-viet@quadram.ac.uk"
LABEL org.opencontainers.image.url='https://github.com/quadram-institute-bioscience/ribosomal_snakemake'
LABEL org.opencontainers.image.documentation='https://github.com/quadram-institute-bioscience/ribosomal_snakemake/readme.md'
LABEL org.opencontainers.image.title='RiboTree'
LABEL org.opencontainers.image.description='Reconstructing phylogenetic tree of bacterial genomes from ribosomal proteins'

SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get -qq -y install wget bzip2

RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc

COPY workflow /opt/workflow

ENV PATH /opt:/opt/conda/bin:/opt/conda/envs/rc2_snakemake/bin:$PATH

RUN source ~/.bashrc \
    && conda install -c conda-forge mamba \
    && mamba env create -f /opt/workflow/envs/environment.yaml \
    && mamba clean --all 

# install tiny init for docker entrypoint
RUN apt-get install -y curl grep sed dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

COPY ribotree.py /opt/

COPY entrypoint.sh /usr/local/bin/entrypoint.sh

RUN chmod 0775 /usr/local/bin/entrypoint.sh

WORKDIR /data

ENV SNAKEMAKE_OUTPUT_CACHE=/data
ENV HOME=/data

ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]

CMD ["/bin/bash"]