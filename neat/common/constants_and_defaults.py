import numpy as np
from frozendict import frozendict

LOW_COVERAGE_THRESHOLD = 50
LOW_PROBABILITY_THRESHOLD = 1e-12

"""
Constants needed for analysis
"""
# allowed nucleotides sets not only the allowed letters, but their order for sampling purposes.
ALLOWED_NUCL = ['A', 'C', 'G', 'T']
NUC_IND = frozendict({'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3})

# IUPAC ambiguity codes mapped to the concrete bases they represent. NEAT resolves each
# occurrence to one of these bases (chosen with the run's seeded RNG) when the reference is
# loaded, so the rest of the pipeline — reads, BAM, golden VCF, error/trinucleotide lookups —
# only ever sees A/C/G/T/N and never crashes on an unmapped base. 'N' is intentionally NOT
# listed: it keeps its existing handling (masked to a low-quality telomere-repeat filler in
# Read.convert_masking) rather than being resolved to a definite base.
IUPAC_CODES = frozendict({
    'R': ('A', 'G'), 'Y': ('C', 'T'), 'S': ('C', 'G'), 'W': ('A', 'T'),
    'K': ('G', 'T'), 'M': ('A', 'C'),
    'B': ('C', 'G', 'T'), 'D': ('A', 'G', 'T'), 'H': ('A', 'C', 'T'), 'V': ('A', 'C', 'G'),
})


def resolve_iupac_bases(sequence: str, rng) -> tuple[str, int]:
    """
    Replace IUPAC ambiguity codes in an (uppercased) sequence with concrete bases.

    Each ambiguity code (R, Y, S, W, K, M, B, D, H, V) is replaced with one of the bases it
    represents, chosen uniformly with the run's seeded RNG. Standard bases (A/C/G/T) and 'N'
    are left untouched — 'N' retains its downstream low-quality masking. Reference assemblies
    such as GRCh38 carry a handful of these codes; resolving them at load time guarantees the
    rest of the pipeline only sees A/C/G/T/N and never raises on an unmapped base.

    The common case (no ambiguity codes present) returns the input unchanged after a single
    vectorized scan, so genome-scale references are not penalized.

    :param sequence: An uppercased nucleotide sequence.
    :param rng: The run's seeded numpy Generator.
    :return: (resolved_sequence, number_of_bases_resolved)
    """
    arr = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    code_bytes = np.fromiter((ord(c) for c in IUPAC_CODES), dtype=np.uint8)
    mask = np.isin(arr, code_bytes)
    if not mask.any():
        return sequence, 0

    arr = arr.copy()
    n_resolved = 0
    for code, bases in IUPAC_CODES.items():
        positions = np.where(arr == ord(code))[0]
        if positions.size == 0:
            continue
        base_bytes = np.fromiter((ord(b) for b in bases), dtype=np.uint8)
        arr[positions] = base_bytes[rng.integers(len(bases), size=positions.size)]
        n_resolved += int(positions.size)
    return arr.tobytes().decode("ascii"), n_resolved

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
