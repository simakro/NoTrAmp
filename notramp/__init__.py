#Copyright (c) 2022, Simon Magin (simakro)
#BSD-2 license
#see LICENSE file for license details

#__doc__=
"""
NoTrAmp is a Tool for read-depth normalization and trimming of amplicon
 reads generated with long read technologies (ONT/PacBio).
"""

# import os
from .version import __version__
import map_trim
import amp_cov
