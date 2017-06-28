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
from skimage import io
import sys
sys.path.append('/n/denic_lab/Users/nweir/python_packages/')
sys.path.append(
    '/n/denic_lab/Users/nweir/python_packages/csth-imaging/dependencies/')
from csth_analysis import czi_io, find_cells, segment_cells, foci

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
# load .czi file into MultiFinder instance
finder = find_cells.MultiFinder(czi_path)
# load bg file from multi-image .czi and add to finder
bg_im_path = '/n/denic_lab/Lab/TH_Imaging/6.8.17 TH69 LSM Airyscan 880/HEK LC3 WIPI for Nick/LC3-WIPI_HEK_dSQSTM1_NoTreat_ClumpandEmptyandSingles_Airyscan Processing.czi'
bg_czi = czi_io.load_multi_czi(bg_im_path)
bg_czi_im = np.expand_dims(bg_czi[0][4, :, :, :, :], axis=0)
finder.bg_im = bg_czi_im
finder.bg_channels = bg_czi[1]
# initialize a CellSplitter from finder
splitter = segment_cells.CellSplitter(finder)
splitter.segment_nuclei(verbose=True)  # segment nuclei
splitter.segment_cells(488, verbose=True)  # segment cells using the 488 wl
# initialize a Foci instance from splitter
foci_obj = foci.Foci(splitter, verbose=True)
foci_obj.segment(verbose=True)  # segment foci using PexSegmenter
foci_obj.count_foci(verbose=True)  # count foci
foci_obj.measure_overlap(verbose=True)  # measure # of overlapping foci
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
foci_obj.pandas_output(output_dir + '/' + str(array_no) + '.csv',
                       verbose=True)
# output images to check quality of segmentation later
im_fname = foci_obj.filenames.split('/')[-1]
im_output_dir = output_dir + '/' + im_fname[:-4]
if not os.path.isdir(im_output_dir):
    os.makedirs(im_output_dir)
os.chdir(im_output_dir)
for i in range(0, len(foci_obj.segmented_nuclei)):
    io.imsave(str(i)+'_nuclei.tif', foci_obj.segmented_nuclei[i])
    io.imsave(str(i)+'_cells.tif', foci_obj.segmented_cells[i])
    for c in foci_obj.foci.keys():  # keys are channel ints
        io.imsave(str(i)+'_'+str(c)+'_foci.tif', foci_obj.foci[c][i])
for c in finder.cell_channels:
    ch_ims = finder.get_channel_arrays(c, bg=False)
    for i in range(0, ch_ims.shape[0]):
        io.imsave(str(i)+'_'+str(c)+'_raw.tif', ch_ims[i, :, :, :])
