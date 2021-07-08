#!/bin/bash --login
set -euo pipefail


conda activate rc2_snakemake

# /usr/bin/tini --

exec "$@"
