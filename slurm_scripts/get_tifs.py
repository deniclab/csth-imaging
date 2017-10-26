import os
import sys
import numpy as np
import argparse
from skimage import io
sys.path.append('/n/denic_lab/Users/nweir/python_packages/')
from csth_analysis import czi_io



parser = argparse.ArgumentParser(description='Extract .tifs from czi files.')
parser.add_argument('-d', '--czi_dir', required = True,
                    help = 'directory containing images to segment.')
parser.add_argument('-o', '--out_dir', required = True,
                    help = 'directory to output tif files to.')

args = parser.parse_args()
print(args)
czi_dir = args.czi_dir
if not czi_dir.endswith('/'):
    czi_dir = czi_dir + '/'
output_dir = args.out_dir
os.chdir(czi_dir)
if not output_dir.endswith('/'):
    output_dir = output_dir + '/'
czi_list = [f for f in os.listdir() if '.czi' in f]
for f in czi_list:
    if 'Empty' in f:
        continue
    print('current file: ' + f)
    czi_arr, channel_arr = czi_io.load_multi_czi(czi_dir + f)
    print('czi array shape:')
    print(czi_arr.shape)
    if len(czi_arr.shape) == 4:
        czi_arr = np.expand_dims(czi_arr, axis=0)
    for im in range(0, czi_arr.shape[0]):
        for c in range(0, czi_arr.shape[1]):
            if channel_arr[c] == 561:
                channel = '594'
            else:
                channel = str(channel_arr[c])
            io.imsave(output_dir+f[:-4]+'_'+str(im)+'_'+channel+'.tif',
                      czi_arr[im, c, :, :, :])
