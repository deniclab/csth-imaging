#!/usr/bin/env
# -*- coding: utf-8 -*-
"""Classes and methods for finding cells in images based on background."""

import czifile
import numpy as np
from warnings import warn
from skimage import io, measure
from csth_analysis import czi_io
from scipy import stats
from scipy.ndimage import filters
from warnings import warn


class CellMask:
    """Container for cell mask and intermediates."""

    # NOTE: I DON'T USE THIS FOR ANYTHING RIGHT NOW!
    def __init__(self, raw_im, raw_bg, gaussian_im,
                 pvals_im, cell_labs, cell_mask):
        """Create an instance of a cell mask.

        NOTE: This class is not intended to be called directly, but should
        instead be initialized by the CellFinder.find_cells() method.
        """
        self.raw_im = raw_im
        self.raw_bg = raw_bg
        self.gaussian_im = gaussian_im
        self.pvals_im = pvals_im
        self.cell_labs = cell_labs
        self.cell_mask = cell_mask


class MultiFinder:
    """Distinguish cells from background in multi-position czi files."""

    def __init__(self, filename, bg_index=-1, bg_filename=''):
        """Create a MultiFinder object.

        Arguments:
            filename (str): The path to an image file to initialize the object.
            bg_index (int, optional): If the initializing file is .czi format,
                this indicates the index within the czi file array that
                corresponds to the background image with no cells. Defaults to
                -1, which means the background image isn't contained within
                the initializing image file.
            bg_filename (str, optional): The path to an image file
                corresponding to the background image. This argument must be
                provided if the background stage position is not contained
                in the multi-image file provided by filename or an
                exception will occur. It should not be provided if the
                background image is contained within a multi-image file
                provided by filename.
        """
        self.filenames = [filename]  # use only fname, not path
        if bg_index == -1:
            if bg_filename == '':
                warn('No background image provided during initialization.')
            self.bg_origin = 'separate'  # separate czi or tiff file
            self.bg_filename = bg_filename
        else:
            self.bg_origin = 'slice'  # slice of a multi-czi
        if '.tif' in self.filenames[0]:
            self.cell_im = io.imread(self.filenames[0])
            # TODO: implement adding dimensions for multi-img/channels both
            # here and a method to add them later
        elif '.czi' in self.filenames[0]:
            cell_czi = czi_io.load_multi_czi(self.filenames[0])
            self.cell_im = cell_czi[0]
            # if cell_im has shape C-Z-X-Y, not IMG-C-Z-X-Y, add axis for img
            if len(self.cell_im.shape) == 4:
                self.cell_im = np.expand_dims(self.cell_im, axis=0)
            # add filename:slice #s dict to indicate which imgs came from where
            self.f_to_s = {self.filenames[0]: range(0, self.cell_im.shape[0])}
            self.cell_channels = cell_czi[1]
        if self.bg_filename != '':
            if '.tif' in self.bg_filename:
                self.bg_im = io.imread(self.bg_filename)
                # TODO: implement adding dimensions for multi-img/channels both
                # here and a method to add them later
            elif '.czi' in self.bg_filename:
                bg_czi = czi_io.load_single_czi(self.bg_filename)
                self.bg_im = np.expand_dims(bg_czi[0], axis=0)
                self.bg_channels = bg_czi[1]
        elif bg_index != -1:
            self.bg_im = self.cell_im[bg_index, :, :, :, :]
            # remove parts of cell_im that correspond to bg
            bg_mask = np.ones(shape=self.cell_im.shape, dtype=bool)
            bg_mask[bg_index, :, :, :, :] = False
            self.cell_im = self.cell_im[bg_mask]
            self.bg_channels = self.cell_channels

    def add_czi(self, filename):
        """Add an additional czi file containing additional image(s).

        Arguments:
            filename (str): Path to the czi file to be added to the existing
                MultiFinder object.
        """
        new_czi = czi_io.load_multi_czi(filename)
        stripped_fname = filename.split('/')[1]
        if len(new_czi[0].shape) == 4:  # if it's a single img czi, not multi
            new_czi[0] = np.expand_dims(new_czi[0], axis=0)
        if new_czi[0].shape[1:] != self.cell_im.shape[1:]:  # if diff CZXY
            raise ValueError('The new czi has a different shape than cell_im.')
        if new_czi[1] != self.cell_channels:
            raise ValueError('The new czi uses non-matching channels.')
        self.filenames.append(stripped_fname)
        self.f_to_s[stripped_fname] = range(
            self.cell_im.shape[0], self.cell_im.shape[0] + new_czi.shape[0])
        self.cell_im = np.concatenate(self.cell_im, new_czi)

    def find_cells(self, channel, return_all=False, verbose=True):
        """Find cells within all images in the indicated channel."""
        # get channel images first
        im_arrs = self.get_channel_arrays(channel)
        # transform into log space, as bg is roughly log-normal
        # this requires adding 1 to each value to avoid NaN log-xform
        if verbose:
            print('log-transforming arrays...')
        log_f_im = np.log10(im_arrs[0] + 1)
        log_bg_im = np.log10(im_arrs[1] + 1)
        if verbose:
            print('applying gaussian filter...')
        log_gaussian_f = filters.gaussian_filter(log_f_im,
                                                 sigma=[0, 0, 3, 3])
        log_gaussian_bg = filters.gaussian_filter(log_bg_im,
                                                  sigma=[0, 0, 3, 3])
        bg_mean = np.mean(log_gaussian_bg)
        bg_sd = np.std(log_gaussian_bg)
        # get p-val that px intensity could be brighter than the value in each
        # array position in the "positive" im, which will indicate where
        # fluorescence is.
        if verbose:
            print('computing p-value transformation...')
        f_pvals = np.empty_like(log_gaussian_f)
        for s in range(0, f_pvals.shape[0]):
            print(' computing p-val xform for slice ' + str(s + 1) +
                  ' out of ' + str(f_pvals.shape[0]))
            f_pvals[s, :, :, :] = 1-stats.norm.cdf(
                log_gaussian_f[s, :, :, :], bg_mean, bg_sd)
        # convert to binary using empirically tested cutoffs (p<0.5/65535)
        f_pvals = f_pvals*65535
        f_pvals = f_pvals.astype('uint16')
        f_pvals_binary = np.copy(f_pvals)
        if verbose:
            print('converting to binary...')
        f_pvals_binary[f_pvals > 0] = 0
        f_pvals_binary[f_pvals == 0] = 1
        # eliminate too-small regions that don't correspond to cells
        cell_masks = []
        if return_all:
            raw_labs_list = []
            labs_list = []
        if verbose:
            print('')
            print('generating cell masks...')
            print('')
        for im in range(0, im_arrs[0].shape[0]):
            if verbose:
                print('generating mask #' + str(im + 1))
            curr_im = f_pvals_binary[im, :, :, :]
            if verbose:
                print('labeling contiguous objects...')
            r_labs = measure.label(curr_im, connectivity=2,
                                   background=0)
            # next command eliminates objects w vol < 100,000 px and
            # generates binary array indicating where cells are for output
            if verbose:
                print('eliminating objects w/volume < 100,000 px...')
            # don't count zeros in next line to avoid including background
            objs_w_cts = np.unique(r_labs[r_labs != 0], return_counts=True)
            cell_mask = np.reshape(np.in1d(
                r_labs, objs_w_cts[0][objs_w_cts[1] > 100000]),
                                   r_labs.shape)
            if verbose:
                print('pruning labels...')
            trim_labs = np.copy(r_labs)
            trim_labs[np.invert(cell_mask)] = 0  # eliminate small obj labels
            if verbose:
                print('appending outputs...')
            cell_masks.append(cell_mask)
            if return_all:
                raw_labs_list.append(r_labs)
                labs_list.append(trim_labs)
            if verbose:
                print('mask #' + str(im + 1) + ' complete.')
                print()
        if return_all:
            return({'input_ims': im_arrs[0],
                    'input_bg': im_arrs[1],
                    'gaussian_ims': np.power(10, log_gaussian_f),
                    'pvals_ims': f_pvals,
                    'raw_labs': raw_labs_list,
                    'trimmed_labs': labs_list,
                    'cell_masks': cell_masks})
        else:
            return cell_masks
    # helper methods #

    def get_channel_arrays(self, channel, fluorescence=True, bg=True,
                           mode='multi', ind=0):
        """Extract im arrays for specific channels."""
        channel = int(channel)  # Implicitly checks that channel is an int
        return_vals = []
        # return tuple of (fluorescence_array, bg_array) for the channel
        if fluorescence:
            if mode == 'multi':
                return_vals.append(self.cell_im[
                    :, self.cell_channels.index(channel), :, :, :])
            elif mode == 'single':
                return_vals.append(self.cell_im[
                    ind, self.cell_channels.index(channel), :, :, :])
        if bg:
            return_vals.append(self.bg_im[
                :, self.bg_channels.index(channel), :, :, :])
        if len(return_vals) == 1:
            return return_vals[0]
        else:
            return tuple(return_vals)


class CellFinder:
    """Distinguish cells from background in fluorescence images."""

    def __init__(self, im_filename, bg_im_filename):
        """Create a CellFinder object."""
        # set attributes
        self.filename = im_filename
        self.bg_filename = bg_im_filename
        if '.tif' in self.filename:
            self.cell_im = io.imread(self.filename)
        elif '.czi' in self.filename:
            cell_czi = czi_io.load_single_czi(self.filename)
            self.cell_im = cell_czi[0]
            self.cell_channels = cell_czi[1]
        if '.tif' in self.bg_filename:
            self.bg_im = io.imread(self.bg_filename)
        elif '.czi' in self.bg_filename:
            bg_czi = czi_io.load_single_czi(self.bg_filename)
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
        log_gaussian_f = filters.gaussian_filter(log_f_im, sigma=[0, 2, 2])
        log_gaussian_bg = filters.gaussian_filter(log_bg_im, sigma=[0, 2, 2])
        bg_mean = np.mean(log_gaussian_bg)
        bg_sd = np.std(log_gaussian_bg)
        # get p-val that px intensity could be brighter than the value in each
        # array position in the "positive" im, which will indicate where
        # fluorescence is.
        f_pvals = 1-stats.norm.cdf(log_gaussian_f, bg_mean, bg_sd)
        # convert to binary using empirically tested cutoffs (p<0.5/65535)
        f_pvals = f_pvals*65535
        f_pvals = f_pvals.astype('uint16')
        f_pvals_binary = np.copy(f_pvals)
        f_pvals_binary[f_pvals > 0] = 0
        f_pvals_binary[f_pvals == 0] = 1
        # eliminate too-small regions that don't correspond to cells
        r_labs = measure.label(f_pvals_binary, connectivity=2,
                               background=0)
        # next command eliminates objects w vol < 100,000 px
        trim_labs = np.reshape(np.in1d(
            r_labs, np.unique(r_labs)[np.unique(
                r_labs, return_counts=True)[1] < 100000]))
        # generate binary array indicating where cells are for output
        cell_mask = np.zeros(shape=f_pvals_binary.shape)
        cell_mask[trim_labs != 0] = 1

        return CellMask(im_arrs[0], im_arrs[1], np.power(10, log_gaussian_f),
                        f_pvals, trim_labs, cell_mask)

    # helper methods #

    def get_channel_arrays(self, channel, fluorescence=True, bg=True):
        """Extract im arrays for specific channels."""
        channel = int(channel)  # Implicitly checks that channel is an int
        return_vals = []
        # return tuple of (fluorescence_array, bg_array) for the channel
        if fluorescence:
            return_vals.append(self.cell_im[
                self.cell_channels.index(channel), :, :, :])
        if bg:
            return_vals.append(self.bg_im[
                self.bg_channels.index(channel), :, :, :])
        return tuple(return_vals)
