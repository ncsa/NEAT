"""
Constants used in the NEAT utilities
"""

__all__ = [
    'REF_WHITELIST',
    'VALID_NUCL',
    'VALID_TRINUC',
    'VCF_DEFAULT_POP_FREQ'
]

#why 30?
REF_WHITELIST = [str(n) for n in range(1, 30)] + ['x', 'y', 'X', 'Y', 'mt', 'Mt', 'MT']
REF_WHITELIST += ['chr' + n for n in REF_WHITELIST]

VALID_NUCL = ['A', 'C', 'G', 'T']

VALID_TRINUC = [VALID_NUCL[i] + VALID_NUCL[j] + VALID_NUCL[k] for i in range(len(VALID_NUCL)) for j in
                range(len(VALID_NUCL)) for k in range(len(VALID_NUCL))]

# if parsing a dbsnp vcf, and no CAF= is found in info tag, use this as default val for population freq
VCF_DEFAULT_POP_FREQ = 0.00001