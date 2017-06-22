#!/usr/bin/env
# -*- coding: utf-8 -*-
"""Classes and methods for watershed segmentation of cells using nuclei."""

import numpy as np
import czi_io
import find_cells
from pyto_segmenter import PexSegment

class Nuclei:
    """Container for segmented nuclei."""

# FOR THIS, WHAT DO I NEED TO DO?
# - load images (load_multi_czi)
# - find cells (MultiFinder.find_cells)
# - find nuclei (PexSegment?)
# - generate distance map for cells to nuclei
# - watershed segment cells
