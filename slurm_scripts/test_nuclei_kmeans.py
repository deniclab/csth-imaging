import os
from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
import argparse
import sys
sys.path.append('/n/denic_lab/Users/nweir/python_packages/')
sys.path.append(
    '/n/denic_lab/Users/nweir/python_packages/csth-imaging/dependencies/')
from csth_analysis import find_cells

parser = argparse.ArgumentParser(description='Process LC3/WIPI stain imgs.')
parser.add_argument('-r', '--ref_csv', required=True,
                    help='path to reference csv file')
parser.add_argument('-a', '--array_no', required=True,
                    help='SLURM job array ID')
parser.add_argument('-o', '--output_dir', required=True,
                    help='dir for CSV-formatted output')
args = parser.parse_args()
print('args:')
print(args)
ref_csv = args.ref_csv
array_no = int(args.array_no)
output_dir = args.output_dir
# read .czi file path from csv reference table
ref_df = pd.read_csv(ref_csv)
czi_path = ref_df['files'].iloc[array_no]
print('czi path: ' + czi_path)
finder = find_cells.MultiFinder(czi_path)
nuclei_ims = finder.get_channel_arrays(405, bg=False)
