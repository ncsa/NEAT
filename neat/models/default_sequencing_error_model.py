"""
The following is the default sequencing error, based on Illumina analysis
"""

import numpy as np
from ..variants import Insertion, Deletion, SingleNucleotideVariant

# The following default model parameters are here, for tweaking and revision.
default_avg_seq_error = 0.009228843915252066
default_read_length = 151

# This puts a high probability toward getting a maximum quality score. The current values
# should be considered temporary. We're working on final values.
chance_of_change = 0.5
chance_of_lower = 0.5
chance_of_higher = 0.5

default_error_variant_probs = {Insertion: 0.004, Deletion: 0.006, SingleNucleotideVariant: 0.99}
default_indel_len_model = {1: 0.999, 2: 0.001}
default_insertion_model = np.array([0.25, 0.25, 0.25, 0.25])
