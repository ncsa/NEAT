"""
The following model is the default for cancer mutation in humans.
"""

import numpy as np

from .default_mutation_model import default_mutation_sub_matrix, default_trinuc_trans_matrices, \
    default_trinuc_trans_bias, default_insertion_lengths, default_deletion_lengths, default_deletion_weights

# The following CANCER default model parameters are here, for tweaking and revision.
# Note that these are not yet implemented
default_cancer_avg_mut_rate = 0.002
default_cancer_homozygous_freq = 0.2
default_cancer_insertion_chance = 0.03
default_cancer_deletion_chance = 0.07
default_cancer_mutation_sub_matrix = default_mutation_sub_matrix
default_cancer_trinuc_trans_matrices = default_trinuc_trans_matrices
default_cancer_trinuc_trans_bias = default_trinuc_trans_bias
default_cancer_insertion_lengths = default_insertion_lengths
default_cancer_insertion_weights = np.array([0.1, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
default_cancer_deletion_lengths = default_deletion_lengths
default_cancer_deletion_weights = default_deletion_weights
