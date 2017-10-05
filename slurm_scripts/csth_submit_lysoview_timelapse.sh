#!/bin/bash

#SBATCH -n 1
#SBATCH -t 0-24:00
#SBATCH -p serial_requeue
#SBATCH --mem=50000
#SBATCH -o /n/denic_lab/Users/nweir/python_packages/csth-imaging/logs/%A_%a.out
#SBATCH -e /n/denic_lab/Users/nweir/python_packages/csth-imaging/logs/%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nweir@fas.harvard.edu

r_csv=$1
out_dir=$2
source new-modules.sh
source activate PYTO_SEG_ENV

python3 /n/denic_lab/Users/nweir/python_packages/csth-imaging/slurm_scripts/csth_submit_lc3_lysoview.py -r $r_csv \
    -a $SLURM_ARRAY_TASK_ID -o $out_dir
