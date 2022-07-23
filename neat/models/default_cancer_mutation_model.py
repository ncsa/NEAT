"""
The following model is the default for cancer mutation in humans.
"""

from .default_mutation_model import default_mutation_sub_matrix, default_trinuc_trans_matrices, \
    default_trinuc_trans_bias, default_deletion_len_model
from ..variants import Insertion, Deletion, SingleNucleotideVariant

# The following CANCER default model parameters are here, for tweaking and revision.
# Note that these are not yet implemented
default_cancer_avg_mut_rate = 0.002
default_cancer_homozygous_freq = 0.2
# See the variants package for orders of types and definitions.
default_cancer_variant_probs = {Insertion: 0.03, Deletion: 0.07, SingleNucleotideVariant: 0.9}

default_cancer_insert_len_model = {1: 0.125, 2: 0.125, 3: 0.3125, 4: 0.0625,
                                   6: 0.0625, 7: 0.0625, 8: 0.0625, 9: 0.0625,
                                   10: 0.0625}
