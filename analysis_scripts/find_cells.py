#!/usr/bin/env
# -*- coding: utf-8 -*-
"""Class and methods for finding cells in images based on background."""

import czifile
import numpy as np
from warnings import warn
from skimage import io
import czi_io
from scip import stats

class CellMask:
    """Container for cell mask and intermediates."""

    pass


class CellFinder:
    """Distinguish cells from background in fluorescence images."""

    def __init__(self, im_filename, bg_im_filename):
        # set attributes
        self.filename = im_filename
        self.bg_filename = bg_im_filename
        if '.tif' in self.filename:
            self.cell_im = io.imread(self.filename)
        elif '.czi' in self.filename:
            cell_czi = czi_io.load_czi(self.filename)
            self.cell_im = cell_czi[0]
            self.cell_channels = cell_czi[1]
        if '.tif' in self.bg_filename:
            self.bg_im = io.imread(self.bg_filename)
        elif '.czi' in self.bg_filename:
            bg_czi = czi_io.load_czi(self.bg_filename)
            self.bg_im = bg_czi[0]
            self.bg_channels = bg_czi[0]
        # check inputs
        if self.filename.shape[0] != self.bg_filename.shape[0]:
            warn('bg image and cell images have different #s of channels.')

    def find_cells(self, channel):
        """Find cells within image in the indicated channel."""

        # get channel images first
        im_arrs = self.get_channel_arrays(channel)
        # transform into log space, as bg is roughly log-normal
        log_f_im = np.log10(im_arrs[0])
        log_bg_im = np.log10(im_arrs[1])
        bg_mean = np.mean(log_bg_im)
        bg_sd = np.std(log_bg_im)
        # get p-val that px intensity could be brighter than the value in each
        # array position in the "positive" im, which will indicate where
        # fluorescence is.
        f_pvals = 1-stats.norm.cdf(log_f_im, bg_mean, bg_sd)
        # TODO NEXT: find sizes of contiguous regions, remove small ones



    # helper methods #

    def get_channel_arrays(self, channel, fluorescence=True, bg=True):
        """A helper method for extracting im arrays for specific channels."""

        channel = int(channel)  # Implicitly checks that channel is an int
        return_vals = []
        # return tuple of (fluorescence_array, bg_array) for the channel
        if fluorescence:
            return_vals.append(self.cell_im[np.argwhere(self.cell_channels ==
                                                        channel), :, :, :])
        if bg:
            return_vals.append(self.bg_im[np.argwhere(self.bg_channels ==
                                                      channel), :, :, :])
        return tuple(return_vals)
