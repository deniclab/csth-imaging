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
    im_path, bg_filename=bg_path, log_path=output_dir + '/log', foc_channel=640)
print('MultiFinder created.')
# initialize a CellSplitter from finder
splitter = segment_cells.CellSplitter(finder, lo_p=True, threshold=300)
print('CellSplitter instance created.')
splitter.segment_nuclei(verbose=True)  # segment nuclei
print('Nuclei segmented.')
splitter.segment_cells(488, verbose=True)  # segment cells using the 488 wl
print('Cells segmented.')
# initialize a Foci instance from splitter
foci_obj = foci.Foci(splitter, verbose=True)
print('Foci instance created.')
foci_obj.segment(seg_channels=(640,), min_cutoff={640: 5}, rm_nuclear=False,
                 thresholds={640: (4000, 2000)})
print('Foci segmented.')
cell_im = foci_obj.segmented_cells[0]
all_foci_mask = foci_obj.foci[640][0] > 0
cell_golgi_vals = pd.DataFrame(index=np.unique(cell_im),
                               columns=['intensity_mean', 'volume'])
for c in np.unique(cell_im):
    if c == 0:
        continue
    cell_foci_mask = np.logical_and(all_foci_mask, cell_im == c)
    golgi_volume = np.sum(cell_foci_mask)
    golgi_mean = np.sum(foci_obj.imgs[488][0][cell_foci_mask])/np.sum(cell_foci_mask)
    cell_golgi_vals.at[c, 'intensity_mean'] = golgi_mean
    cell_golgi_vals.at[c, 'volume'] = golgi_volume
cell_golgi_vals.to_csv(output_dir +h 'golgi.csv')
# output images to check quality of segmentation later
print('outputting images...')
im_fname = foci_obj.filenames.split('/')[-1]
im_output_dir = output_dir + '/' + im_fname[:-4]
if not os.path.isdir(im_output_dir):
    os.makedirs(im_output_dir)
os.chdir(im_output_dir)
for i in range(0, len(foci_obj.segmented_nuclei)):
    skimage.io.imsave(str(i)+'_nuclei.tif', foci_obj.segmented_nuclei[i].astype('uint16'))
    skimage.io.imsave(str(i)+'_cells.tif', foci_obj.segmented_cells[i].astype('uint16'))
    for c in foci_obj.foci.keys():  # keys are channel ints
        skimage.io.imsave(str(i)+'_'+str(c)+'_foci.tif', foci_obj.foci[c][i].astype('uint16'))
for c in finder.cell_channels:
    ch_ims = finder.get_channel_arrays(c, bg=False)
    for i in range(0, ch_ims.shape[0]):
        skimage.io.imsave(str(i)+'_'+str(c)+'_raw.tif', ch_ims[i, :, :, :].astype('uint16'))
