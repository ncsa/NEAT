"""
The following is the default sequencing error, based on Illumina analysis
"""

import numpy as np
from ..variants import Insertion, Deletion, SingleNucleotideVariant

# The following default model parameters are here, for tweaking and revision.
default_avg_seq_error = 0.01
default_read_length = 151

# This is based on the input model NEAT was using previously
default_error_transition_matrix = np.array(
    [[0.0, 0.4918, 0.3377, 0.1705],
     [0.5238, 0.0, 0.2661, 0.2101],
     [0.3755, 0.2355, 0.0, 0.389],
     [0.2505, 0.2552, 0.4943, 0.0]]
)

# This list may not be the final list
default_quality_scores = np.array([2, 11, 25, 37])

# This puts a high probability toward getting a maximum quality score. The current values
# should be considered temporary. We're working on final values.
default_qual_score_probs = np.array()

default_error_variant_probs = {Insertion: 0.004, Deletion: 0.006, SingleNucleotideVariant: 0.99}
default_indel_len_model = {1: 0.999, 2: 0.001}
default_insertion_model = np.array([0.25, 0.25, 0.25, 0.25])
