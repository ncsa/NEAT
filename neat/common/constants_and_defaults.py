from frozendict import frozendict

LOW_COVERAGE_THRESHOLD = 50
LOW_PROBABILITY_THRESHOLD = 1e-12

"""
Constants needed for analysis
"""
# allowed nucleotides sets not only the allowed letters, but their order for sampling purposes.
ALLOWED_NUCL = ['A', 'C', 'G', 'T']
NUC_IND = frozendict({'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3})

MAX_ATTEMPTS = 100  # max attempts to insert a mutation into a valid position
MAX_MUTFRAC = 0.3  # the maximum percentage of a window that can contain mutations

ALL_TRINUCS = [ALLOWED_NUCL[i] + ALLOWED_NUCL[j] + ALLOWED_NUCL[k]
               for i in range(len(ALLOWED_NUCL))
               for j in range(len(ALLOWED_NUCL))
               for k in range(len(ALLOWED_NUCL))]
ALL_CONTEXTS = [f'{ALLOWED_NUCL[i]}_{ALLOWED_NUCL[k]}'
                for i in range(len(ALLOWED_NUCL))
                for k in range(len(ALLOWED_NUCL))]
TRINUC_IND = frozendict({ALL_TRINUCS[i]: i for i in range(len(ALL_TRINUCS))})
DINUC_IND = frozendict({ALL_CONTEXTS[i]: i for i in range(len(ALL_CONTEXTS))})

# DEBUG
IGNORE_TRINUC = False

# percentile resolution used for fraglen quantizing
COV_FRAGLEN_PERCENTILE = 10.
LARGE_NUMBER = 9999999999

# For humans. Whitelist is used to generate the mutation model
HUMAN_WHITELIST = [str(n) for n in range(1, 30)] + ['x', 'y', 'X', 'Y', 'mt', 'Mt', 'MT']
HUMAN_WHITELIST += ['chr' + n for n in HUMAN_WHITELIST]