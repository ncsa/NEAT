import copy
import pathlib

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

ALLOWED_NUCL = ['A', 'C', 'G', 'T']
OK_CHR_ORD = {'A': True, 'C': True, 'G': True, 'T': True, 'U': True}
TRI_IND = {'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3, 'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
           'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11, 'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15}
NUC_IND = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
ALL_TRI = [ALLOWED_NUCL[i] + ALLOWED_NUCL[j] + ALLOWED_NUCL[k] for i in range(len(ALLOWED_NUCL)) for j in range(len(ALLOWED_NUCL)) for k in range(len(ALLOWED_NUCL))]
ALL_IND = {ALL_TRI[i]: i for i in range(len(ALL_TRI))}

# DEBUG
IGNORE_TRINUC = False

# percentile resolution used for fraglen quantizing
COV_FRAGLEN_PERCENTILE = 10.
LARGE_NUMBER = 9999999999

"""
DEFAULT MUTATION MODELS
TODO: Convert to pickle files? Why are these done different?

Convert to namespace?

class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

model_1 = Namespace(
                    DEFAULT_1_OVERALL_MUT_RATE = 0.001,
                    DEFAULT_1_HOMOZYGOUS_FREQ = 0.010,
                    DEFAULT_1_INDEL_FRACTION = 0.05,
                    DEFAULT_1_INS_VS_DEL = 0.6,
                    DEFAULT_1_INS_LENGTH_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    DEFAULT_1_INS_LENGTH_WEIGHTS = [0.4, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.034, 0.033, 0.033],
                    DEFAULT_1_DEL_LENGTH_VALUES = [1, 2, 3, 4, 5],
                    DEFAULT_1_DEL_LENGTH_WEIGHTS = [0.3, 0.2, 0.2, 0.2, 0.1],
                    example_matrix_1 = [[0.0, 0.15, 0.7, 0.15],
                                        [0.15, 0.0, 0.15, 0.7],
                                        [0.7, 0.15, 0.0, 0.15],
                                        [0.15, 0.7, 0.15, 0.0]],
                    DEFAULT_1_TRI_FREQS = [copy.deepcopy(example_matrix_1) for _ in range(16)],
                    DEFAULT_1_TRINUC_BIAS = [1. / float(len(ALL_TRI)) for _ in ALL_TRI],
                    DEFAULT_MODEL_1 = [DEFAULT_1_OVERALL_MUT_RATE,
                                       DEFAULT_1_HOMOZYGOUS_FREQ,
                                       DEFAULT_1_INDEL_FRACTION,
                                       DEFAULT_1_INS_VS_DEL,
                                       DEFAULT_1_INS_LENGTH_VALUES,
                                       DEFAULT_1_INS_LENGTH_WEIGHTS,
                                       DEFAULT_1_DEL_LENGTH_VALUES,
                                       DEFAULT_1_DEL_LENGTH_WEIGHTS,
                                       DEFAULT_1_TRI_FREQS,
                                       DEFAULT_1_TRINUC_BIAS]
                    )
                    
model_2 = Namespace(
                    DEFAULT_2_OVERALL_MUT_RATE = 0.002,
                    DEFAULT_2_HOMOZYGOUS_FREQ = 0.200,
                    DEFAULT_2_INDEL_FRACTION = 0.1,
                    DEFAULT_2_INS_VS_DEL = 0.3,
                    DEFAULT_2_INS_LENGTH_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    DEFAULT_2_INS_LENGTH_WEIGHTS = [0.1, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05],
                    DEFAULT_2_DEL_LENGTH_VALUES = [1, 2, 3, 4, 5],
                    DEFAULT_2_DEL_LENGTH_WEIGHTS = [0.3, 0.2, 0.2, 0.2, 0.1],
                    example_matrix_2 = [[0.0, 0.15, 0.7, 0.15],
                                        [0.15, 0.0, 0.15, 0.7],
                                        [0.7, 0.15, 0.0, 0.15],
                                        [0.15, 0.7, 0.15, 0.0]],
                    DEFAULT_2_TRI_FREQS = [copy.deepcopy(example_matrix_2) for _ in range(16)],
                    DEFAULT_2_TRINUC_BIAS = [1. / float(len(ALL_TRI)) for _ in ALL_TRI],
                    DEFAULT_MODEL_2 = [DEFAULT_2_OVERALL_MUT_RATE,
                                       DEFAULT_2_HOMOZYGOUS_FREQ,
                                       DEFAULT_2_INDEL_FRACTION,
                                       DEFAULT_2_INS_VS_DEL,
                                       DEFAULT_2_INS_LENGTH_VALUES,
                                       DEFAULT_2_INS_LENGTH_WEIGHTS,
                                       DEFAULT_2_DEL_LENGTH_VALUES,
                                       DEFAULT_2_DEL_LENGTH_WEIGHTS,
                                       DEFAULT_2_TRI_FREQS,
                                       DEFAULT_2_TRINUC_BIAS]
                    )
"""

DEFAULT_1_OVERALL_MUT_RATE = 0.001
DEFAULT_1_HOMOZYGOUS_FREQ = 0.010
DEFAULT_1_INDEL_FRACTION = 0.05
DEFAULT_1_INS_VS_DEL = 0.6
DEFAULT_1_INS_LENGTH_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
DEFAULT_1_INS_LENGTH_WEIGHTS = [0.4, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.034, 0.033, 0.033]
DEFAULT_1_DEL_LENGTH_VALUES = [1, 2, 3, 4, 5]
DEFAULT_1_DEL_LENGTH_WEIGHTS = [0.3, 0.2, 0.2, 0.2, 0.1]
example_matrix_1 = [[0.0, 0.15, 0.7, 0.15],
                    [0.15, 0.0, 0.15, 0.7],
                    [0.7, 0.15, 0.0, 0.15],
                    [0.15, 0.7, 0.15, 0.0]]
DEFAULT_1_TRI_FREQS = [copy.deepcopy(example_matrix_1) for _ in range(16)]
DEFAULT_1_TRINUC_BIAS = [1. / float(len(ALL_TRI)) for _ in ALL_TRI]
DEFAULT_MODEL_1 = [DEFAULT_1_OVERALL_MUT_RATE,
                   DEFAULT_1_HOMOZYGOUS_FREQ,
                   DEFAULT_1_INDEL_FRACTION,
                   DEFAULT_1_INS_VS_DEL,
                   DEFAULT_1_INS_LENGTH_VALUES,
                   DEFAULT_1_INS_LENGTH_WEIGHTS,
                   DEFAULT_1_DEL_LENGTH_VALUES,
                   DEFAULT_1_DEL_LENGTH_WEIGHTS,
                   DEFAULT_1_TRI_FREQS,
                   DEFAULT_1_TRINUC_BIAS]

DEFAULT_2_OVERALL_MUT_RATE = 0.002
DEFAULT_2_HOMOZYGOUS_FREQ = 0.200
DEFAULT_2_INDEL_FRACTION = 0.1
DEFAULT_2_INS_VS_DEL = 0.3
DEFAULT_2_INS_LENGTH_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
DEFAULT_2_INS_LENGTH_WEIGHTS = [0.1, 0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
DEFAULT_2_DEL_LENGTH_VALUES = [1, 2, 3, 4, 5]
DEFAULT_2_DEL_LENGTH_WEIGHTS = [0.3, 0.2, 0.2, 0.2, 0.1]
example_matrix_2 = [[0.0, 0.15, 0.7, 0.15],
                    [0.15, 0.0, 0.15, 0.7],
                    [0.7, 0.15, 0.0, 0.15],
                    [0.15, 0.7, 0.15, 0.0]]
DEFAULT_2_TRI_FREQS = [copy.deepcopy(example_matrix_2) for _ in range(16)]
DEFAULT_2_TRINUC_BIAS = [1. / float(len(ALL_TRI)) for _ in ALL_TRI]
DEFAULT_MODEL_2 = [DEFAULT_2_OVERALL_MUT_RATE,
                   DEFAULT_2_HOMOZYGOUS_FREQ,
                   DEFAULT_2_INDEL_FRACTION,
                   DEFAULT_2_INS_VS_DEL,
                   DEFAULT_2_INS_LENGTH_VALUES,
                   DEFAULT_2_INS_LENGTH_WEIGHTS,
                   DEFAULT_2_DEL_LENGTH_VALUES,
                   DEFAULT_2_DEL_LENGTH_WEIGHTS,
                   DEFAULT_2_TRI_FREQS,
                   DEFAULT_2_TRINUC_BIAS]
