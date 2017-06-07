#!/usr/bin/env
# -*- coding: utf-8 -*-
"""Class and methods for finding cells in images based on background."""

import czifile
import numpy as np
from skimage import io
import czi_io


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
            self.cell_im = czi_io.load_czi(self.filename)[0]
        if '.tif' in self.bg_filename:
            self.bg_im = io.imread(self.bg_filename)
        elif '.czi' in self.bg_filename:
            self.bg_im = czi_io.load_czi(self.bg_filename)
        # check inputs
        if self.filename.shape != self.bg_filename.shape:
            raise ValueError('bg image and cell images are different shapes.')
