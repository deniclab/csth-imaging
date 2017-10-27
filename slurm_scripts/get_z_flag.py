import os
import numpy as np
import pandas as pd
import argparse
from skimage import io
import pickle
import sys
sys.path.append('/n/denic_lab/Users/nweir/python_packages/')
sys.path.append(
    '/n/denic_lab/Users/nweir/python_packages/csth-imaging/dependencies/')
from csth_analysis import find_cells

parser = argparse.ArgumentParser(description='Process LC3/p62 stain imgs.')
parser.add_argument('-r', '--ref_csv', required=True,
                    help='path to reference csv file')
parser.add_argument('-o', '--output_path', required=True,
                    help='path for CSV-formatted output')
args = parser.parse_args()
print('args:')
print(args)
ref_csv = args.ref_csv
array_no = int(args.array_no)
output_path = args.output_path
# read .czi file path from csv reference table
ref_df = pd.read_csv(ref_csv)
czi_paths = ref_df['files'].tolist()
# load .czi file into MultiFinder instance
output_df = pd.DataFrame(columns=['filename', 'im_number', 'flagged_z'])
for czi in czi_paths:
    finder = find_cells.MultiFinder(
        czi,
        oof_svm='/n/denic_lab/Users/nweir/python_packages/csth-imaging/trained_svm.pkl')
    im_for_clf = finder.get_channel_arrays(finder.foc_channel, bg=False)
    print('unlabeling out of focus slices...')
    with open(finder.oof_svm, 'rb') as r:
        clf = pickle.load(r)
    shrt_fname = finder.filenames[0].split('/')[-1][:-4]
    for im in range(0, im_for_clf.shape[0]):
        if finder.log_path is not None:
            focus_slices = finder.get_blur_slices(
                im=im_for_clf[im, :, :, :], clf=clf, slc_no=im,
                log_path=finder.log_path+'/'+shrt_fname)
        else:
            focus_slices = finder.get_blur_slices(
                im=im_for_clf[im, :, :, :], clf=clf, slc_no=im)
        if focus_slices[0] == 1 or focus_slices[-1] == 1:
            finder.flagged_z_ims[im] = 1
    curr_df = pd.DataFrame({'filename': finder.filenames,
                            'image': list(range(0, finder.cell_im.shape[0])),
                            'flagged_z': finder.flagged_z_ims
                            })
    output_df = pd.concat([output_df, curr_df], ignore_index=True)
output_df.to_csv(output_path, index=False)
