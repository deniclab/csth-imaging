#!/bin/bash

#SBATCH -n 1
#SBATCH -t 0-3:00
#SBATCH -p serial_requeue
#SBATCH --mem=20000
#SBATCH -o /n/denic_lab/Users/nweir/python_packages/csth-imaging/logs/%j.out # STDOUT
#SBATCH -e /n/denic_lab/Users/nweir/python_packages/csth-imaging/logs/%j.err # STDERR
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nweir@fas.harvard.edu

src_dir=$1
dest_dir=$2

source new-modules.sh
source activate PYTO_SEG_ENV

python3 /n/denic_lab/Users/nweir/python_packages/csth-imaging/slurm_scripts/get_tifs.py -d $src_dir -o $dest_dir
