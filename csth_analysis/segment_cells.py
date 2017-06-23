#!/usr/bin/env
# -*- coding: utf-8 -*-
"""Classes and methods for watershed segmentation of cells using nuclei."""

import numpy as np
import czi_io
import find_cells
from pyto_segmenter import PexSegment
from scipy.ndimage.morphology import distance_transform_edt
from scipy.ndimage.filters import gaussian_filter
from skimage.morphology import watershed


class CellSplitter:
    """Class and methods for segmenting cells using watershedding from DAPI."""

    def __init__(self, multi_finder, high_threshold=1500,
                 low_threshold=750):
        """Create a Nuclei object for segmentation."""
        self.multi_finder = multi_finder
        self.high_threshold = high_threshold
        self.low_threshold = low_threshold
        self.segmented_nuclei = []
        self.nuclei_centers = []
        if 405 not in multi_finder.cell_channels:
            raise ValueError(
                'The MultiFinder object lacks nuclei fluorescence.'
            )
        self.seg_input = self.multi_finder.get_channel_arrays(405)[0]
        if len(self.seg_input.shape) == 4:
            self.nimgs = self.seg_input.shape[0]
            self.input_shape = 4
        else:
            if len(self.seg_input.shape) == 3:
                self.input_shape = 3
            else:
                raise ValueError('The 405 wl image is an unexpected shape.')

    def segment_nuclei(self):
        """Segment nuclei from 405 nm wavelength images."""
        if self.input_shape == 3:
            segmenter = PexSegment(
                src_data=self.seg_input, mode='canny',
                high_threshold=self.high_threshold,
                low_threshold=self.low_threshold)
            seg_output = segmenter.segment()
            self.segmented_nuclei.append(seg_output.peroxisomes)
            self.nuclei_centers.append(seg_output.maxima)
        elif self.input_shape == 4:
            for i in range(0, self.seg_input.shape[0]):
                segmenter = PexSegment(
                    src_data=self.seg_input[i, :, :, :], mode='canny',
                    high_threshold=self.high_threshold,
                    low_threshold=self.low_threshold)
                seg_output = segmenter.segment()
                self.segmented_nuclei.append(seg_output.peroxisomes)
                self.nuclei_centers.append(seg_output.maxima)
        # remove perinuclear foci that are smaller than a true nucleus
        for i in range(0, len(self.segmented_nuclei)):
            objs_w_cts = np.unique(self.segmented_nuclei[i],
                                   return_counts=True)
            objs_to_rm = objs_w_cts[0][objs_w_cts[1] < 1000]
            self.segmented_nuclei[i][np.reshape(
                np.in1d(self.segmented_nuclei[i], objs_to_rm),
                self.segmented_nuclei[i].shape)] = 0
            # remove corresponding nuclei_centers as well
            self.nuclei_centers[i][self.segmented_nuclei[i] == 0] = 0

    def segment_cells(self, channel):
        """Segment cells, identified using find_cells, based on nuclei."""
        segmented_cells = []
        # test to make sure nuclei have already been segmented.
        if len(self.segmented_nuclei) == 0:
            self.segment_nuclei()
        self.cell_masks = self.multi_finder.find_cells(channel, verbose=False)
        # convert segmented nuclei to an inverted mask for distance xform
        nuclei_masks = self.segmented_nuclei
        for i in nuclei_masks:
            i[i > 0] = 1
            i = np.invert(i.astype('bool'))
        # segment cells
        for j in range(0, len(self.cell_masks)):
            # distance xform based on the distance to a nucleus
            dist_map = distance_transform_edt(nuclei_masks[j],
                                              sampling=(3, 1, 1))
            dist_map = gaussian_filter(dist_map, [1, 2, 2])  # smooth the map
            # generate segmentation seeds from nuclei segmentation maxima
            labs = PexSegment.PexSegmenter.watershed_labels(
                self.nuclei_centers[j])
            # watershed segment
            segmented_cells.append(
                watershed(dist_map, labs, mask=self.cell_masks[j]))


# FOR THIS, WHAT DO I NEED TO DO?
# - load images (load_multi_czi)
# - find cells (MultiFinder.find_cells)
# - find nuclei (PexSegment?)
# - generate distance map for cells to nuclei
# - watershed segment cells
