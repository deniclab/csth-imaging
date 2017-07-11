from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
import argparse
import sys
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
sys.path.append('/n/denic_lab/Users/nweir/python_packages/')
sys.path.append(
    '/n/denic_lab/Users/nweir/python_packages/csth-imaging/dependencies/')
from csth_analysis import find_cells

parser = argparse.ArgumentParser(description='Process LC3/WIPI stain imgs.')
parser.add_argument('-r', '--ref_csv', required=True,
                    help='path to reference csv file')
parser.add_argument('-a', '--array_no', required=True,
                    help='SLURM job array ID')
args = parser.parse_args()
print('args:')
print(args)
ref_csv = args.ref_csv
array_no = int(args.array_no)
# read .czi file path from csv reference table
ref_df = pd.read_csv(ref_csv)
czi_path = ref_df['files'].iloc[array_no]
print('czi path: ' + czi_path)
finder = find_cells.MultiFinder(czi_path)
nuclei_ims = finder.get_channel_arrays(405, bg=False)
print('nuclei_ims shape:')
print(nuclei_ims.shape)
del finder  # save memory
cm = plt.cm.get_cmap('RdYlBu_r')  # colormap for plotting
for i in range(0, nuclei_ims.shape[0]):
    c_im = nuclei_ims[i, :, :, :].astype('uint16')
    kmeans_subset = np.random.choice(c_im.flatten(), 1000000)
    cluster_out = KMeans(n_clusters=2).fit(kmeans_subset.reshape(-1, 1))
    n, bins = np.histogram(kmeans_subset, bins=100)
    colors_1 = np.zeros(shape=bins.shape)
    colors_1[bins > np.amin(
        kmeans_subset.reshape(-1, 1)[cluster_out.labels_ == 1])] = 1
    colors_1 = colors_1[:-1]
    plt.bar(bins[:-1], n, color=cm(colors_1), width=bins[1]-bins[0])
    plt.yscale('log', nonposy='clip')
    fname = czi_path.split('/')[-1]
    fname = fname[:-4] + '_' + str(i) + '.pdf'
    plt.savefig('/n/denic_lab/Lab/csth-output/nuclei_kmeans_test/' + fname)
