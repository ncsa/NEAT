from frozendict import frozendict

LOW_COVERAGE_THRESHOLD = 50
LOW_PROBABILITY_THRESHOLD = 1e-12

"""
Constants needed for analysis
"""
# allowed nucleotides sets not only the allowed letters, but their order for sampling purposes.
ALLOWED_NUCL = ['A', 'C', 'G', 'T']
NUC_IND = frozendict({'A': 0, 'C': 1, 'G': 2, 'T': 3})

MAX_ATTEMPTS = 100  # max attempts to insert a mutation into a valid position
MAX_MUTFRAC = 0.3  # the maximum percentage of a window that can contain mutations

DINUC_IND = frozendict({'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3, 'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
                        'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11, 'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15})
ALL_TRI = [ALLOWED_NUCL[i] + ALLOWED_NUCL[j] + ALLOWED_NUCL[k] for i in range(len(ALLOWED_NUCL))
           for j in range(len(ALLOWED_NUCL)) for k in range(len(ALLOWED_NUCL))]
TRI_IND = frozendict({ALL_TRI[i]: i for i in range(len(ALL_TRI))})

# DEBUG
IGNORE_TRINUC = False

# percentile resolution used for fraglen quantizing
COV_FRAGLEN_PERCENTILE = 10.
LARGE_NUMBER = 9999999999

# For humans. Whitelist is used to generate the mutation model
HUMAN_WHITELIST = [str(n) for n in range(1, 30)] + ['x', 'y', 'X', 'Y', 'mt', 'Mt', 'MT']