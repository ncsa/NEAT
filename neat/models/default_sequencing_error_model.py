"""
The following is the default sequencing error, based on Illumina analysis
"""

import numpy as np

# The following default model parameters are here, for tweaking and revision.
default_avg_seq_error = 0.01
default_read_length = 151

# This is based on the input model NEAT was using previously
default_transition_matrix = np.array(
    [[0.0, 0.4918, 0.3377, 0.1705],
     [0.5238, 0.0, 0.2661, 0.2101],
     [0.3754, 0.2355, 0.0, 0.389],
     [0.2505, 0.2552, 0.4942, 0.0]]
)

# This list may not be the final list
default_quality_scores = np.array([0, 12, 24, 36])

# This puts a high probability toward getting a maximum quality score. The current values
# should be considered temporary. We're working on final values.
default_qual_score_probs = np.array([
    np.full(151, 0.00001),
    np.full(151, 0.001),
    np.full(151, 0.23),
    np.full(151, 0.76899)
]).T

default_insertion_probability = 0.004
default_deletion_probability = 0.006
default_indel_lengths = np.array([1, 2])
default_indel_weights = np.array([0.999, 0.001])
default_insertion_model = np.array([0.25, 0.25, 0.25, 0.25])
