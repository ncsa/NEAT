"""
The following model is the default human mutation model.
"""

import numpy as np

from ..common import ALL_TRI

# The following default model parameters are here, for tweaking and revision.
default_avg_mut_rate = 0.001
default_homozygous_freq = 0.001
default_insertion_chance = 0.03
default_deletion_chance = 0.02
default_mutation_sub_matrix = np.array(
    [[0.0, 0.15, 0.7, 0.15],
     [0.15, 0.0, 0.15, 0.7],
     [0.7, 0.15, 0.0, 0.15],
     [0.15, 0.7, 0.15, 0.0]]
)
default_trinuc_trans_matrices = np.full((16, 4, 4), default_mutation_sub_matrix)
default_trinuc_trans_bias = np.full(len(ALL_TRI), 1 / float(len(ALL_TRI)))
default_insertion_lengths = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
default_insertion_weights = np.array([0.4, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.034, 0.033, 0.033])
default_deletion_lengths = np.array([1, 2, 3, 4, 5])
default_deletion_weights = np.array([0.3, 0.2, 0.2, 0.2, 0.1])
