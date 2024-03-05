"""
Constants used in the NEAT utilities
"""

import numpy as np


__all__ = [
    'VCF_DEFAULT_POP_FREQ',
    'DEF_HOMOZYGOUS_FRQ',
    'DEF_MUT_SUB_MATRIX'
]

# if parsing a dbsnp vcf, and no CAF= is found in info tag, use this as default val for population freq
VCF_DEFAULT_POP_FREQ = 0.00001

#make homozygous constant and trans matrix
DEF_HOMOZYGOUS_FRQ = 0.001

DEF_MUT_SUB_MATRIX = np.array(
    [[0.0, 0.15, 0.7, 0.15],
     [0.15, 0.0, 0.15, 0.7],
     [0.7, 0.15, 0.0, 0.15],
     [0.15, 0.7, 0.15, 0.0]]
)

#DEF_TRINUC_TRANS_MATRIX = np.full((16, 4, 4), DEF_MUT_SUB_MATRIX)
