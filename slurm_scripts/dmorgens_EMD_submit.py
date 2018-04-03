# steps for this script:
# 1. get arguments:
#   - path to .csv ref file
#   - array ID
#   - path to output directory
#
# 2. use pd.read_csv to read ref file
# 3. open .czi according to path
# 4. perform analysis
# 5. save output csv to output dir
import os
import numpy as np
import pandas as pd
import argparse
import skimage
import sys
sys.path.append('/n/denic_lab/Users/nweir/python_packages/')
sys.path.append(
    '/n/denic_lab/Users/nweir/python_packages/csth-imaging/dependencies/')
from csth_analysis import find_cells, segment_cells, foci

parser = argparse.ArgumentParser(description='Process LC3/p62 stain imgs.')
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
im_path = ref_df['files'].iloc[array_no]
bg_path = ref_df['controls'].iloc[array_no]
print('im path: ' + im_path)
print('control path: ' + bg_path)
# load .czi file into MultiFinder instance
finder = find_cells.MultiFinder(
    im_path, bg_filename=bg_path, log_path=output_dir + '/log', foc_channel=488)
print('MultiFinder created.')
# initialize a CellSplitter from finder
splitter = segment_cells.CellSplitter(finder, threshold=200)
print('CellSplitter instance created.')
splitter.segment_nuclei(verbose=True)  # segment nuclei
print('Nuclei segmented.')
dilation_structure = np.array([[[0,0,0],
                                [0,0,0],
                                [0,0,0]],
                               [[0,1,0],
                                [1,1,1],
                                [0,1,0]],
                               [[0,0,0],
                                [0,0,0],
                                [0,0,0]]])
dilated_labeled_boundaries = np.zeros(shape=splitter.segmented_nuclei[0].shape)
print('Unique boundaries: {}'.format(np.unique(splitter.segmented_nuclei[0])))
for n in np.unique(splitter.segmented_nuclei[0]):
    print('Processing boundary {}'.format(n))
    if n == 0:
        continue
    current_nucleus = splitter.segmented_nuclei[0] == n
    current_nucleus = ndi.morphology.binary_erosion(
        current_nucleus, iterations=10, structure=dilation_structure)
    current_boundary = np.zeros(shape=current_nucleus.shape)
    for s in range(0, current_nucleus.shape[0]):
        current_boundary[s, :, :] = skimage.segmentation.find_boundaries(
            current_nucleus[s, :, :], mode='inner')
    dilated_boundary = ndi.morphology.binary_dilation(
        current_boundary, iterations=4, structure=dilation_structure)
    dilated_labeled_boundaries[dilated_boundary == 1] = n
skimage.io.imsave(
    output_dir+'/log/'+splitter.filenames[0][:-4]+'_boundaries.tif',
    dilated_labeled_boundaries)
nuc_ids = []
nuc_means = []
print('Calculating nuclear boundary fluorescences.')
for n in np.unique(dilated_labeled_boundaries):
    if n == 0:
        continue
    mean_by_slice = np.empty(shape=dilated_labeled_boundaries.shape[0])
    for s in range(mean_by_slice.size):
        mean_by_slice[s] = np.mean(
            splitter.multi_finder.cell_im[0, 0, s, :, :][
                dilated_labeled_boundaries[s] == n])
    mean_by_slice = np.sort(mean_by_slice[~np.isnan(mean_by_slice)])[::-1]
    nuc_border_mean = np.mean(mean_by_slice[0:5])
    print('Nucleus {} mean: {}'.format(n, nuc_border_mean))
    nuc_ids.append(n)
    nuc_means.append(nuc_border_mean)
print('Outputting data to csv.')
output_df = pd.DataFrame({'nucleus_id': nuc_ids,
                          'intensity_mean': nuc_border_mean})
output_df.to_csv(output_dir + splitter.filenames[0][:-4]+'.csv')
