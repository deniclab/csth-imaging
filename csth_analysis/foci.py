#!/usr/bin/env
# -*- coding: utf-8 -*-
"""Classes and methods for segmenting and assigning intracellular foci."""

from pyto_segmenter.PexSegment import PexSegmenter
import numpy as np
import pandas as pd
from scipy.ndimage.morphology import binary_erosion


class Foci:
    """Container for segmented foci."""

    def __init__(self, CellSplitter, verbose=True):
        """Create a Foci instance from a CellSplitter instance."""
        if verbose:
            print('initializing attributes...')
        self.filenames = CellSplitter.filenames
        if len(self.filenames) == 1:
            self.filenames = self.filenames[0].split('/')[-1]
        self.segmented_nuclei = CellSplitter.segmented_nuclei
        self.segmented_cells = CellSplitter.segmented_cells
        self.n_raw_nuclei = CellSplitter.n_raw_nuclei
        self.cell_masks = CellSplitter.cell_masks
        self.n_cells = CellSplitter.n_cells
        if verbose:
            print('loading images...')
        self.imgs = {}
        for c in CellSplitter.multi_finder.cell_channels:
            self.imgs[c] = CellSplitter.multi_finder.get_channel_arrays(
                c, bg=False)
        self.channels = CellSplitter.multi_finder.cell_channels
        self.n_pos = self.imgs[self.channels[0]].shape[0]  # num of stage posns

    def segment(self, verbose=True):
        """Identify foci in image."""
        self.foci = {}
        if verbose:
            print('beginning segmentation.')
        self.thresholds = {488: 12500, 561: 6000}  # diff thresh for 488/561
        self.erosion_struct = np.array(  # the strel for erosion of nuclei
            [[[False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, True,  False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False]],
             [[False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, True,  False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False]],
             [[False, False, False,  True, False, False, False],
              [False, False,  True,  True,  True, False, False],
              [False,  True,  True,  True,  True,  True, False],
              [True,   True,  True,  True,  True,  True,  True],
              [False,  True,  True,  True,  True,  True, False],
              [False, False,  True,  True,  True, False, False],
              [False, False, False,  True, False, False, False]],
             [[False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, True,  False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False]],
             [[False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, True,  False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False],
              [False, False, False, False, False, False, False]]])
        for c in self.channels:
            if c == 405:
                continue  # don't segment foci from DAPI channel
            channel_foci = []
            if verbose:
                print('------------------------------------------------------')
                print('segmenting foci from channel ' + str(c))
                print('------------------------------------------------------')
                print('canny threshold for channel ' + str(c) + ': ' +
                      str(self.thresholds[c]))
            for i in range(0, self.n_pos):  # for each stage position
                # segment foci from this channel
                if verbose:
                    print('segmenting foci for position ' + str(i + 1) +
                          ' out of ' + str(self.n_pos))
                curr_segmenter = PexSegmenter(
                    src_data=self.imgs[c][i, :, :, :], seg_method='canny',
                    high_threshold=self.thresholds[c],
                    low_threshold=self.thresholds[c]/2)
                curr_seg = curr_segmenter.segment()
                c_foci = curr_seg.peroxisomes
                raw_img = curr_seg.raw_img
                # get #s and volumes for foci and make dict
                if verbose:
                    print('eliminating dim foci...')
                objs, vols = np.unique(c_foci, return_counts=True)
                if verbose:
                    print('before eliminating dim foci: ' +
                          str(len(objs) - 1) + ' foci in image')
                vols = dict(zip(objs, vols))
                # get mean intensity for each focus
                mean_intensity = {}
                for obj in objs:
                    mean_intensity[obj] = np.sum(
                        raw_img[c_foci == obj]).astype('float')/vols[obj]
                rev_dict = {v: k for k, v in mean_intensity.items()}
                cell_mean = np.mean(raw_img[self.segmented_cells[i] != 0])
                cell_sd = np.std(raw_img[self.segmented_cells[i] != 0])
                for k in rev_dict.keys():
                    if k < cell_mean + 3*cell_sd:  # if mean intensity too low
                        c_foci[c_foci == rev_dict[k]] = 0  # eliminate focus
                if verbose:
                    print('after eliminating dim foci: ' +
                          str(len(np.unique(c_foci))-1) + ' foci in image')
                if verbose:
                    print('eliminating foci that reside outside of cells...')
                c_foci[self.cell_masks[i] == 0] = 0
                if verbose:
                    print('eliminating intranuclear foci...')
                eroded_nuclei = np.copy(self.segmented_nuclei[i])
                eroded_nuclei = binary_erosion(eroded_nuclei,
                                               structure=self.erosion_struct)
                c_foci[eroded_nuclei != 0] = 0
                if verbose:
                    print(str(len(np.unique(c_foci))-1) + ' final foci')
                channel_foci.append(c_foci)
                if verbose:
                    print('foci segmented from position ' + str(i + 1))
                    print()
            self.foci[c] = channel_foci
            if verbose:
                print('------------------------------------------------------')
                print('segmentation complete for channel ' + str(c))
                print('------------------------------------------------------')

    def count_foci(self, verbose=True):
        """Count the # of foci present in each segmented cell.

        Yields:
            a dict of dicts.
            Inner dict: a dict of image # (key): ndarray of #s of foci per cell
                (val) pairs.
            Outer dict: a dict of channel (key): Inner dict (val) pairs.
                channel is represented as a string to allow later addition of
                'overlap' as a key through self.measure_overlap().

        """
        if not hasattr(self, 'foci'):
            if verbose:
                print('foci not yet segmented. calling segment().')
            self.segment()
        self.foci_cts = {}
        for k in self.foci.keys():  # for each channel
            if verbose:
                print('------------------------------------------------------')
                print('counting foci/cell in channel ' + str(k) + ' images...')
                print('------------------------------------------------------')
            im_foci = {}
            for i in range(0, len(self.foci[k])):  # for each image
                print('counting foci/cell in image #' + str(i+1) + ' of ' +
                      str(len(self.foci[k])))
                foci_per_cell = []
                for cell in np.unique(self.segmented_cells[i]):
                    if cell == 0:
                        continue  # skip background
                    # get unique segmented foci IDs within the cell
                    cell_foci = np.unique(
                        self.foci[k][i][self.segmented_cells[i] == cell])
                    foci_per_cell.append(len(cell_foci) - 1)  # subtract bg
                im_foci[i] = np.asarray(foci_per_cell)
            self.foci_cts[str(k)] = im_foci

    def measure_overlap(self, channel=1, verbose=True):
        """Determine how many foci overlap between two channels."""
        if not hasattr(self, 'foci_cts'):
            if verbose:
                print('foci have not been counted yet. calling count_foci().')
            self.count_foci()
        if verbose:
            print('checking channels...')
        if verbose:
            print('checking for overlap...')
        channels = list(self.foci.keys())
        if channel == 1:
            overlap_channel = channels[0]
        else:
            overlap_channel = channels[1]
        if verbose:
            print('measuring overlap in channel ' + str(overlap_channel))
        n_ims = len(self.foci[overlap_channel])
        overlap = {}
        self.overlapping_foci = []
        self.overlap_channel = overlap_channel
        im_overlaps = {}
        for i in range(0, n_ims):
            if verbose:
                print('finding overlap in image #' + str(i + 1) + 'out of ' +
                      str(n_ims))
            overlap = np.logical_and(self.foci[channels[0]][i] > 0,
                                     self.foci[channels[1]][i] > 0)
            if verbose:
                print('getting IDs of overlapping foci...')
            overlap_IDs = np.unique(self.foci[overlap_channel][i][overlap])
            overlap_foci = np.copy(self.foci[overlap_channel][i])
            if verbose:
                print('creating image of overlapping foci...')
            overlap_foci[np.reshape(
                np.in1d(overlap_foci, overlap_IDs, invert=True),
                overlap_foci.shape)] = 0
            self.overlapping_foci.append(overlap_foci)
            if verbose:
                print('counting # of overlapping foci per cell...')
            cell_overlap_n = []
            for c in np.unique(self.segmented_cells[i]):
                if c == 0:
                    continue  # skip background
                cell_overlap_n.append(len(
                    np.unique(overlap_foci[self.segmented_cells[i] == c])) - 1)
            im_overlaps[i] = np.asarray(cell_overlap_n)
        self.foci_cts['overlap'] = im_overlaps

    def pandas_output(self, path, verbose=True):
        """Write # of foci and overlap data to a .csv file."""
        if not hasattr(self, 'foci_cts'):
            raise AttributeError('# of foci not yet counted.')
        if 'overlap' not in list(self.foci_cts.keys()):
            raise ValueError('overlap not measured between channels.')
        if verbose:
            print('initializing pd.DataFrame for output...')
        channels = [str(k) for k in self.foci.keys()]
        # initialize arrays that will populate the DataFrame
        tot_foci_1 = np.array([])
        tot_foci_2 = np.array([])
        overlap_foci = np.array([])
        cell_nums = np.array([])
        raw_cells = np.array([])
        total_cells = np.array([])
        im_nums = np.array([])
        if verbose:
            print('populating output arrays...')
        for i in range(0, len(self.segmented_cells)):
            n_cells = len(self.foci_cts[channels[0]][i])
            tot_foci_1 = np.concatenate(
                (tot_foci_1, self.foci_cts[channels[0]][i]))
            tot_foci_2 = np.concatenate(
                (tot_foci_2, self.foci_cts[channels[1]][i]))
            overlap_foci = np.concatenate(
                (overlap_foci, self.foci_cts['overlap'][i]))
            cell_nums = np.concatenate(
                (cell_nums, np.arange(1, n_cells + 1)))
            im_nums = np.concatenate(
                (im_nums, np.repeat(i+1, n_cells)))
            raw_cells = np.concatenate(
                (raw_cells, np.repeat(self.n_raw_nuclei[i], n_cells)))
            total_cells = np.concatenate(
                (total_cells, np.repeat(n_cells, n_cells)))
        non_olap_1 = tot_foci_1 - overlap_foci
        non_olap_2 = tot_foci_2 - overlap_foci
        output_df = pd.DataFrame({'filename': self.filenames,
                                  'image': im_nums,
                                  'cell': cell_nums,
                                  'im_n_cells': total_cells,
                                  'im_raw_cells': raw_cells,
                                  str(channels[0]) + '_total_foci': tot_foci_1,
                                  str(channels[1]) + '_total_foci': tot_foci_2,
                                  str(channels[0]) + '_only_foci': non_olap_1,
                                  str(channels[1]) + '_only_foci': non_olap_2,
                                  'overlap_foci': overlap_foci
                                  })
        output_df.to_csv(path)
