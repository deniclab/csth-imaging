#!/usr/bin/env
# -*- coding: utf-8 -*-
"""Classes and methods for watershed segmentation of cells and nuclei.

Classes:
    CellSplitter: A container for MultiFinder image objects.
        Static methods:
        -rm_edge_objs: eliminate segmented objects that touch the array border.
"""
# Import analysis packages
import numpy as np
from pyto_segmenter.PexSegment import PexSegmenter
import scipy.ndimage.morphology as morph
from scipy.ndimage.morphology import distance_transform_edt
from scipy.ndimage.filters import gaussian_filter
from skimage.morphology import watershed


class CellSplitter(object):
    """
    Class and methods for segmenting cells using watershedding from DAPI.

    Methods:
        __init__: Create a CellSplitter instance from a MultiFinder instance.
        segment_nuclei: Segment nuclei using a 405 channel image.
        segment_cells: Segment cells using nuclei and masks from a MultiFinder.
    Static methods:
        rm_edge_objs: eliminate segmented objects that contact an array border.

    """

    def __init__(self, multi_finder, channel=488, threshold='auto',
                 pval_threshold=0, lo_p=False, cellfinder_mode='pval',
                 cellfinder_threshold=300):
        """
        Create a CellSplitter instance for segmenting cells and nuclei.

        Arguments:
        ----------
        multi_finder : MultiFinder instance
            A MultiFinder created from image files. Contains the images to be
            used for segmentation as well as the methods for generating cell
            masks. Masks don't necessarily need to be generated before creating
            an instance of this class - it will check for them and generate
            them if necessary. See find_cells.py and segment_cells for details.
        channel : int, optional
            Which fluorescence channel to use for segmenting cells. Only used
            if the MultiFinder instance does not already have cell masks made.
            Defaults to 488.
        threshold : {'auto', int}, optional
            How threshold is determined for nuclei segmentation. If 'auto'
            (default), intensity threshold for nuclear segmentation is set
            automatically. If an integer is provided, that value is used as
            the threshold.

        """
        # initialize attributes
        self.filenames = multi_finder.filenames
        self.multi_finder = multi_finder
        self.threshold = threshold
        self.n_raw_nuclei = []
        self.segmented_nuclei = []
        self.nuclei_centers = []
        print('generating cell masks...')
        # generates cell masks using MultiFinder methods from find_cells
        self.cell_masks = self.multi_finder.find_cells(
            channel, verbose=True, pval_threshold=pval_threshold, lo_p=lo_p,
            mode=cellfinder_mode, threshold=cellfinder_threshold)
        if 405 not in multi_finder.cell_channels:  # only takes 405 wl nuclei
            raise ValueError(
                'The MultiFinder object lacks nuclei fluorescence.'
            )
        # extract nuclei images from the MultiFinder
        self.seg_input = self.multi_finder.get_channel_arrays(405)[0]
        # alter shape if only a single image was present in the MultiFinder
        if len(self.seg_input.shape) == 4:
            self.nimgs = self.seg_input.shape[0]
            self.input_shape = 4
        else:
            if len(self.seg_input.shape) == 3:
                self.input_shape = 3
            else:
                raise ValueError('The 405 wl image is an unexpected shape.')

    def segment_nuclei(self, verbose=True):
        """Segment nuclei from 405 nm wavelength images.

        Arguments:
        ---------
        verbose : bool, optional
            Verbose text output. Defaults to True.

        """
        if self.input_shape == 3:  # if only one image to be segmented
            if verbose:
                print('input has 3 dimensions.')
                print('segmenting nuclei using PexSegmenter...')
            if self.threshold == 'auto':
                # use slice-by-slice scaled cutoffs with hard floor to make
                # thresholded image.
                if verbose:
                    print('Using thresholds set based on slice intensity')
                # 1. smooth image
                gaussian_im = gaussian_filter(self.seg_input, sigma=(0, 4, 4))
                # 2. Set the segmentation threshold minimum at 0.2X max px
                thresh_floor = int(np.amax(gaussian_im.flatten())*0.2)
                maxima = np.amax(gaussian_im, axis=0)  # slice-by-slice maxima
                if verbose:
                    print('thresholding floor set to ' + str(thresh_floor))
                    print('slice-by-slice cutoffs:')
                    print(maxima*0.5)  # 50% max slice intensity is threshold
                    print('generating threshold image...')
                threshold_im = np.copy(gaussian_im)
                threshold_im[self.cell_masks == 0] = 0
                for z in range(0, threshold_im.shape[0]):
                    if maxima[z]*0.5 < thresh_floor:
                        # use 0.2x stack max cutoff if slice 0.5x max is low
                        # set all px below the cutoff to 0
                        threshold_im[z, :, :][
                            threshold_im[z, :, :] < thresh_floor] = 0
                    else:  # if slice max*0.5 > stack max*0.2
                        # set all px below cutoff to 0
                        threshold_im[z, :, :][
                            threshold_im[z, :, :] < int(maxima[z]*0.5)] = 0
                threshold_im[threshold_im > 0] = 1  # make bool mask
                # initialize segmenter object for finding nuclei using these
                # masks
                segmenter = PexSegmenter(
                    src_data=threshold_im, seg_method='pre-thresholded')
            else:  # use threshold provided in CellSplitter.__init__
                segmenter = PexSegmenter(
                    src_data=self.seg_input, threshold=self.threshold, g_xy=4,
                    g_z=0)
            # use empirically determined segmentation parameters to segment
            # nuclei
            seg_output = segmenter.segment(fill_holes=True,
                                           edt_sampling=(10, 1, 1),
                                           edt_smooth=[3, 50, 50])
            if verbose:
                print('passing segmented objs and seeds to Nuclei instance...')
            self.segmented_nuclei.append(seg_output.peroxisomes)
            self.nuclei_centers.append(seg_output.maxima)
            del seg_output  # clear up memory
        elif self.input_shape == 4:
            # see pipeline with input_shape == 3 above for clarifying comments
            if verbose:
                print('input image has 4 dimensions (multiple imgs)')
            for i in range(0, self.seg_input.shape[0]):
                if verbose:
                    print('segmenting image #' + str(i+1) +
                          ' of ' + str(self.seg_input.shape[0]) + str('...'))
                if self.threshold == 'auto':
                    # use slice-by-slice scaled cutoffs with hard floor to make
                    # thresholded image.
                    if verbose:
                        print('Using thresholds set based on slice intensity')
                    gaussian_im = gaussian_filter(
                        self.seg_input[i, :, :, :], sigma=(0, 4, 4))
                    thresh_floor = int(np.amax(gaussian_im.flatten())*0.2)
                    # get slice-by-slice max
                    maxima = np.amax(gaussian_im, axis=(1, 2))
                    if verbose:
                        print('thresholding floor set to ' + str(thresh_floor))
                        print('slice-by-slice cutoffs:')
                        print(maxima*0.5)
                        print('generating threshold image...')
                    threshold_im = np.copy(gaussian_im)
                    threshold_im[self.cell_masks[i] == 0] = 0
                    for z in range(0, threshold_im.shape[0]):
                        if maxima[z]*0.5 < thresh_floor:
                            threshold_im[z, :, :][
                                threshold_im[z, :, :] < thresh_floor] = 0
                        else:
                            threshold_im[z, :, :][
                                threshold_im[z, :, :] < int(maxima[z]*0.5)] = 0
                    threshold_im[threshold_im > 0] = 1
                    segmenter = PexSegmenter(
                        src_data=threshold_im, seg_method='pre-thresholded')
                else:
                    if verbose:
                        print('performing segmentation with a threshold of ' +
                              str(self.threshold))
                    segmenter = PexSegmenter(
                        src_data=self.seg_input[i, :, :, :],
                        threshold=self.threshold, g_xy=4, g_z=0)
                seg_output = segmenter.segment(fill_holes=True,
                                               edt_sampling=(10, 1, 1),
                                               edt_smooth=[3, 50, 50])
                if verbose:
                    print('passing segmented objs and seeds' +
                          ' to Nuclei instance...')
                    print()
                self.nuclei_centers.append(seg_output.maxima)
                self.segmented_nuclei.append(seg_output.peroxisomes)
                self.n_raw_nuclei.append(seg_output.npexs)
        # remove perinuclear foci that are smaller than a true nucleus
        if verbose:
            print('removing small perinuclear foci...')
        for i in range(0, len(self.segmented_nuclei)):
            objs_w_cts = np.unique(self.segmented_nuclei[i],
                                   return_counts=True)
            if verbose:
                print(str(len(objs_w_cts[0]) - 1) + ' raw segmented DAPI foci')
                print('object volumes in pixels:')
                print(objs_w_cts[1][1:])
            objs_to_rm = objs_w_cts[0][objs_w_cts[1] < 10000]
            if verbose:
                print('removing objects volume < 10000 px, #s:')
                print(objs_to_rm)
            self.segmented_nuclei[i][np.reshape(
                np.in1d(self.segmented_nuclei[i], objs_to_rm),
                self.segmented_nuclei[i].shape)] = 0
            # remove corresponding nuclei_centers as well
            if verbose:
                print('removing nuclei_centers corresponding to small' +
                      ' DAPI foci...')
            self.nuclei_centers[i][self.segmented_nuclei[i] == 0] = 0

    def segment_cells(self, channel, verbose=True, rm_edge_cells=True):
        """Segment cells, identified using find_cells, based on nuclei.

        Arguments:
        ----------
        channel : int
            Fluorescence channel to use for segmenting cells.
        verbose : bool, optional
            Verbose text output. Defaults to True.
        rm_edge_cells : bool, optional
            Should segmented cells that contact the xy edges of the image be
            removed? Defaults to True.

        """
        strel = np.array([[[0, 0, 0],
                           [0, 0, 0],
                           [0, 0, 0]],
                          [[0, 1, 0],
                           [1, 1, 1],
                           [0, 1, 0]],
                          [[0, 0, 0],
                           [0, 0, 0],
                           [0, 0, 0]]])
        self.segmented_cells = []
        self.n_cells = []
        # test to make sure nuclei have already been segmented.
        if verbose:
            print('checking that nuclei have already been segmented...')
        if len(self.segmented_nuclei) == 0:
            if verbose:
                print('nuclei have not been segmented yet.')
                print('segmenting nuclei...')
            self.segment_nuclei()
        # convert segmented nuclei to an inverted mask for distance xform
        nuclei_masks = np.copy(self.segmented_nuclei)
        channel_ims = self.multi_finder.get_channel_arrays(
            channel, bg=False)
        if verbose:
            print('eroding cell edges...')
        for i in range(0, len(self.cell_masks)):
            if verbose:
                print('eroding mask ' + str(i) + '...')
            means = []
            # use binary erosion to shrink edges of image, finding mean each
            # time. Once erosion begins happening in the cell, the mean won't
            # decrease much more; however, while blurry edge material is being
            # removed by erosion, the cell mean will increase. the inflection
            # point where the slope begins to decrease is called the true
            # edge of the cell, and masks are eroded this many pixels.
            for j in range(0, 10):
                curr_mask = morph.binary_erosion(
                    self.cell_masks[i], structure=strel, iterations=j)
                means.append(
                    np.mean(channel_ims[i, :, :, :][np.logical_and(
                        curr_mask != 0, nuclei_masks[i] == 0)]))
            means[0] = np.mean(channel_ims[i, :, :, :][np.logical_and(
                self.cell_masks[i] != 0, nuclei_masks[i] == 0)])
            if verbose:
                print('cell mask means:')
                print(means)
            slopes = []
            for j in range(0, len(means)-2):
                slopes.append(np.divide(means[j+2]-means[j], 2))
            if verbose:
                print('slopes:')
                print(slopes)
            delta_slopes = np.array([])
            iterations = np.arange(1.5, 8.5, 1)
            for j in range(0, len(slopes)-1):
                delta_slopes = np.append(
                    delta_slopes, np.absolute(slopes[j+1]-slopes[j]))
            delta_slopes[np.isnan(delta_slopes)] = 0
            if verbose:
                print('slope deltas:')
                print(delta_slopes)
            desired_erosions = int(
                iterations[np.argmax(delta_slopes)] + 0.5)
            print('desired erosions: ' + str(desired_erosions))
            self.cell_masks[i] = morph.binary_erosion(
                self.cell_masks[i], structure=strel,
                iterations=desired_erosions)
            if verbose:
                print('erosion of mask #' + str(i) + ' complete.')
        if verbose:
            print('converting nuclei to binary masks for distance xform...')
        for i in range(0, len(nuclei_masks)):
            nuclei_masks[i][nuclei_masks[i] > 0] = 1
            nuclei_masks[i] = np.invert(nuclei_masks[i].astype('bool'))
        # segment cells
        for j in range(0, len(self.cell_masks)):
            if verbose:
                print('segmenting cells in image #' + str(j + 1) +
                      ' out of ' + str(len(self.cell_masks)) + '...')
            # distance xform based on the distance to a nucleus
            if verbose:
                print('performing Euclidean distance xform...')
            dist_map = distance_transform_edt(nuclei_masks[j],
                                              sampling=(3, 1, 1))
            if verbose:
                print('smoothing the distance map...')
            dist_map = gaussian_filter(dist_map, [1, 2, 2])  # smooth the map
            # generate segmentation seeds from nuclei segmentation maxima
            if verbose:
                print('generating segmentation seeds from nuclei_centers...')
            labs = PexSegmenter.watershed_labels(self.nuclei_centers[j])
            # watershed segment
            if verbose:
                print('watershed segmenting cells from nuclei centers...')
                seg_cells = watershed(dist_map, labs, mask=self.cell_masks[j])
            if verbose:
                print('segmentation complete for image #' + str(j + 1) + '.')
            if rm_edge_cells:
                if verbose:
                    print('removing cells that contact the edge of the image.')
                seg_cells = self.rm_edge_objs(seg_cells)
            # identify foci that are smaller than a true cell
            if verbose:
                print('removing small foci...')
            # get tuple of (cell IDs, # of px/cell ID)
            objs_w_cts = np.unique(seg_cells, return_counts=True)
            if verbose:
                print(str(len(objs_w_cts[0]) - 1) + ' raw segmented cells')
            objs_to_rm = objs_w_cts[0][objs_w_cts[1] < 100000]
            print(str(len(objs_to_rm)) +
                  ' foci volume < 100000 px, removing...')
            # remove those objects
            seg_cells[np.reshape(
                np.in1d(seg_cells, objs_to_rm), seg_cells.shape)] = 0
            # calculate stats and pass to class instance
            n = len(np.unique(seg_cells)) - 1  # number of cells
            self.segmented_cells.append(seg_cells)
            self.n_cells.append(n)
            if verbose:
                print()

    @staticmethod
    def rm_edge_objs(arr, z=False):
        """Remove segmented objects that contact the edge of an image.

        Arguments:
        ----------
        arr : 3D np.ndarray
            a 3D NumPy array of ints produced by watershed segmentation.
        z : bool, optional)
            should objects that appear on the top or bottom slice of the
            image be removed (are the first and last slices part of the
            'border'?) defaults to False.

        Returns:
        --------
        An np.ndarray the same shape as `arr` with objects contacting the edge
        changed to 0.

        """
        border_arr = np.ones(shape=arr.shape)
        if z:  # if objects on the top and bottom slices should be removed:
            border_arr[1:-1, 1:-1, 1:-1] = 0  # leave 1s filling top + bottom
        else:
            border_arr[:, 1:-1, 1:-1] = 0  # remove 1s from center of all zs
        border_arr = border_arr.astype('bool')
        # determine which segmented cells contact the edge
        objs_on_edge = np.unique(arr[border_arr])
        output_arr = np.copy(arr)
        # get rid of objects contacting the edge
        output_arr[np.reshape(np.in1d(arr, objs_on_edge), arr.shape)] = 0
        return output_arr
