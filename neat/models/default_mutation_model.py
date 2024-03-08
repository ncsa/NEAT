"""
The following model is the default human mutation model.
"""

import numpy as np

from ..common import ALL_TRINUCS
from ..variants import *

# The following default model parameters are here, for tweaking and revision.
default_avg_mut_rate = 0.001
default_homozygous_freq = 0.001

# See the variants package for orders of types and definitions.
default_variant_probs = {Insertion: 0.03, Deletion: 0.02, SingleNucleotideVariant: 0.95}

default_mutation_sub_matrix = np.array(
    [[0.0, 0.15, 0.7, 0.15],
     [0.15, 0.0, 0.15, 0.7],
     [0.7, 0.15, 0.0, 0.15],
     [0.15, 0.7, 0.15, 0.0]]
)
default_trinuc_trans_matrices = np.full((16, 4, 4), default_mutation_sub_matrix)
default_trinuc_mut_bias = np.full(len(ALL_TRINUCS), 1 / float(len(ALL_TRINUCS)))
default_insertion_len_model = {1: 0.4, 2: 0.2, 3: 0.1, 4: 0.05, 5: 0.05, 6: 0.05, 7: 0.05, 8: 0.034, 9: 0.033, 10: 0.033}
default_deletion_len_model = {1: 0.3, 2: 0.2, 3: 0.2, 4: 0.2, 5: 0.1}
