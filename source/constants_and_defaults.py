import copy

"""
These used to be free-floating in gen_reads. Some of these might be better as inputs. 
"""

"""
Constants needed for runs
"""
VERSION = 4.0
LOW_COVERAGE_THRESHOLD = 50
LOW_PROBABILITY_THRESHOLD = 1e-12

"""
Constants needed for analysis
"""
MAX_ATTEMPTS = 100  # max attempts to insert a mutation into a valid position
MAX_MUTFRAC = 0.3  # the maximum percentage of a window that can contain mutations

# allowed nucleotides sets not only the allowed letters, but their order for sampling purposes.
ALLOWED_NUCL = ['A', 'C', 'G', 'T']
NUC_IND = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

TRI_IND = {'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3, 'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
           'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11, 'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15}
ALL_TRI = [ALLOWED_NUCL[i] + ALLOWED_NUCL[j] + ALLOWED_NUCL[k] for i in range(len(ALLOWED_NUCL)) for j in range(len(ALLOWED_NUCL)) for k in range(len(ALLOWED_NUCL))]
ALL_IND = {ALL_TRI[i]: i for i in range(len(ALL_TRI))}

# DEBUG
IGNORE_TRINUC = False

# percentile resolution used for fraglen quantizing
COV_FRAGLEN_PERCENTILE = 10.
LARGE_NUMBER = 9999999999

# For humans. Whitelist is used to generate the mutation model
HUMAN_WHITELIST = [str(n) for n in range(1, 30)] + ['x', 'y', 'X', 'Y', 'mt', 'Mt', 'MT']

"""
The following are the parameters and models that NEAT uses as a default for its mutation model. Note that I am not
sure where these parameters come from. They were not specific in their naming originally (DEFAULT_MODEL_1), so it
took some investigation to get that these were the mutation model. 
"""
# Not sure where any of these numbers come from
DEFAULT_MUT_MODEL_OVERALL_MUT_RATE = 0.001
DEFAULT_MUT_MODEL_HOMOZYGOUS_FREQ = 0.010
# Percent of mutations that are indels (v snps)
DEFAULT_MUT_MODEL_INDEL_FRACTION = 0.05
# If the mutation is an indel, chance that it is an insertion
DEFAULT_MUT_MODEL_INDEL_INS_PCT = 0.6
# We could do longer insertions simple by making this list longer, but I'm not sure where the weights come from
# They seem to follow some form of 1/x type pattern
DEFAULT_MUT_MODEL_INS_LENGTH_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
DEFAULT_MUT_MODEL_INS_LENGTH_WEIGHTS = [0.4, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.034, 0.033, 0.033]
# Again, longer deletions can be added just by extending this out. Looks like again a 1/x type relationship
DEFAULT_MUT_MODEL_DEL_LENGTH_VALUES = [1, 2, 3, 4, 5]
DEFAULT_MUT_MODEL_DEL_LENGTH_WEIGHTS = [0.3, 0.2, 0.2, 0.2, 0.1]
# This is a subsitution matrix for the set (A,C,T,G). This may be a standard matrix, but I'm not sure where the numbers
# come from but they do have several empirical examples in the paper
DEFAULT_SUBSTITUTION_MATRIX = [[0.0, 0.15, 0.7, 0.15],
                               [0.15, 0.0, 0.15, 0.7],
                               [0.7, 0.15, 0.0, 0.15],
                               [0.15, 0.7, 0.15, 0.0]]

DEFAULT_MUT_MODEL_TRIUC_FREQS = [copy.deepcopy(DEFAULT_SUBSTITUTION_MATRIX) for _ in range(16)]
DEFAULT_MUT_MODEL_TRINUC_BIAS = {x: 1. / float(len(ALL_TRI)) for x in ALL_TRI}


DEFAULT_MUTATION_MODEL = [DEFAULT_MUT_MODEL_OVERALL_MUT_RATE,
                          DEFAULT_MUT_MODEL_HOMOZYGOUS_FREQ,
                          DEFAULT_MUT_MODEL_INDEL_FRACTION,
                          DEFAULT_MUT_MODEL_INDEL_INS_PCT,
                          DEFAULT_MUT_MODEL_INS_LENGTH_VALUES,
                          DEFAULT_MUT_MODEL_INS_LENGTH_WEIGHTS,
                          DEFAULT_MUT_MODEL_DEL_LENGTH_VALUES,
                          DEFAULT_MUT_MODEL_DEL_LENGTH_WEIGHTS,
                          DEFAULT_MUT_MODEL_TRIUC_FREQS,
                          DEFAULT_MUT_MODEL_TRINUC_BIAS]

# This model is used in the cancer sections, so I'm changing the names from DEFAULT_2 to DEFAULT_CANCER_MUTATION_MODEL
DEFAULT_CANCER_MUTATION_OVERALL_MUT_RATE = 0.002
DEFAULT_CANCER_MUTATION_HOMOZYGOUS_FREQ = 0.200
# Fraction of mutations that are indels (instead of SNPs)
DEFAULT_CANCER_MUTATION_INDEL_FRACTION = 0.1
# Percent chance that if the mutation is an indel, that it is an insertion (as opposed to a deletion)
DEFAULT_CANCER_MUTATION_INDEL_INS_PCT = 0.3
DEFAULT_CANCER_MUTATION_INS_LENGTH_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
DEFAULT_CANCER_MUTATION_INS_LENGTH_WEIGHTS = [0.1, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
# Note that all other parameters are identical to default model

DEFAULT_CANCER_MUTATION_MODEL = [DEFAULT_CANCER_MUTATION_OVERALL_MUT_RATE,
                                 DEFAULT_CANCER_MUTATION_HOMOZYGOUS_FREQ,
                                 DEFAULT_CANCER_MUTATION_INDEL_FRACTION,
                                 DEFAULT_CANCER_MUTATION_INDEL_INS_PCT,
                                 DEFAULT_CANCER_MUTATION_INS_LENGTH_VALUES,
                                 DEFAULT_CANCER_MUTATION_INS_LENGTH_WEIGHTS,
                                 DEFAULT_MUT_MODEL_DEL_LENGTH_VALUES,
                                 DEFAULT_MUT_MODEL_DEL_LENGTH_WEIGHTS,
                                 DEFAULT_MUT_MODEL_TRIUC_FREQS,
                                 DEFAULT_MUT_MODEL_TRINUC_BIAS]
