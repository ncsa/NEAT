"""
Creates a model of GC bias in a dataset
"""

import logging
from pathlib import Path

from .utils import get_gc_bias_weights
from ..models import GCBiasModel
from ..common import validate_output_path

__all__ = [
    "compute_gc_bias_runner"
]

_LOG = logging.getLogger(__name__)

def compute_gc_bias_runner(
        bam_file: str | Path,
        reference_file: str | Path,
        output_dir: str | Path,
        output_prefix: str,
        window_size: int = 100,
        overwrite: bool = False):
    """
    Main function for computing GC bias model.
    """
    _LOG.info("Generating GC bias model")

    output_file = Path(output_dir) / (output_prefix + ".pickle.gz")
    validate_output_path(output_file, True, overwrite)

    _LOG.debug(f"BAM: {bam_file}")
    _LOG.debug(f"Reference: {reference_file}")
    _LOG.debug(f"Window size: {window_size}")
    _LOG.debug(f"Output: {output_file}")

    _LOG.info("Calculating GC bias")
    
    # Simple estimation:
    # 1. Calculate GC content for many windows across the reference.
    # 2. Count number of reads mapping to those windows.
    # 3. Normalize to get weights.
    
    weights = get_gc_bias_weights(bam_file, reference_file, window_size)

    model = GCBiasModel(weights, window_size)

    _LOG.info(f"Saving model: {output_file}")
    model.save(output_file)

    _LOG.info("Modeling complete.")
