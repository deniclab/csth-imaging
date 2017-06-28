#!/usr/bin/env
# -*- coding: utf-8 -*-
"""Classes and methods for segmenting and assigning intracellular foci."""

from pyto_segmenter.PexSegment import PexSegmenter
import numpy as np


class Foci:
    """Container for segmented foci."""

    def __init__(self, CellSplitter, verbose=True):
        """Create a Foci instance from a CellSplitter instance."""
        if verbose:
            print('initializing attributes...')
        self.segmented_cells = CellSplitter.segmented_cells
        self.cell_masks = CellSplitter.cell_masks
        self.n_cells = CellSplitter.n_cells
        self.imgs = dict(zip(
            CellSplitter.multi_finder.cell_channels,
            CellSplitter.multi_finder.cell_im))
        self.channels = CellSplitter.multi_finder.cell_channels
        self.n_pos = self.imgs[self.channels[0]].shape[0]  # num of stage posns
