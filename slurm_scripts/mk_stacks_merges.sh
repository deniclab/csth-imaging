#!/bin/bash

#SBATCH -n 1
#SBATCH -t 0-04:00
#SBATCH -p serial_requeue
#SBATCH --mem-per-cpu=4000
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nweir@fas.harvard.edu

img_dir=$1
colors=$2

source new-modules.sh
module load Anaconda3/2.1.0-fasrc01
source activate PYTO_SEG_ENV
module load fiji/1.49j-fasrc01

fiji --headless ~/code/batch_merge_channels.py $img_dir $colors
