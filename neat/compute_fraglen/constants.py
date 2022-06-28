"""
Constants used in the NEAT utilities
"""

__all__ = [
    'FILTER_MINREADS',
    'FILTER_MEDDEV_M',
    'FILTER_MAPQUAL'
]

FILTER_MINREADS = 100  # only consider fragment lengths that have at least this many read pairs supporting it
FILTER_MEDDEV_M = 10  # only consider fragment lengths this many median deviations above the median
FILTER_MAPQUAL = 10  # only consider reads that are mapped with at least this mapping quality
