"""
Constants used in the NEAT utilities
"""

__all__ = [
    'FILTER_MEDDEV_M',
    'FILTER_MAPQUAL'
]

FILTER_MEDDEV_M = 10  # only consider fragment lengths this many median deviations above the median
FILTER_MAPQUAL = 10  # only consider reads that are mapped with at least this mapping quality
