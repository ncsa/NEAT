"""
Emit simulation_summary.json next to the other simulator outputs.

Consumed by `neat compare-vcfs` (issue #297) to attribute false negatives from a
downstream variant caller against the simulator's configuration (mutation bed,
target bed, simulated contigs).
"""
import json
import logging
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pysam

from ... import __version__ as NEAT_VERSION

_LOG = logging.getLogger(__name__)

SCHEMA_VERSION = "1"


def write_simulation_summary(
    options,
    output_dir: Path,
    file_prefix: str,
    config_path: Path,
    analysis_start: float,
    contigs_simulated: Iterable[str],
    reference_contigs: Iterable[str] | None = None,
) -> Path:
    """
    Write `simulation_summary.json` into `output_dir`.

    Returns the path to the written file.
    """
    started_at = _iso_utc(analysis_start)
    completed_ts = time.time()
    completed_at = _iso_utc(completed_ts)

    vcf_path = _abs_or_none(getattr(options, "vcf", None))
    bam_path = _abs_or_none(getattr(options, "bam", None))
    fq1_path = _abs_or_none(getattr(options, "fq1", None))
    fq2_path = _abs_or_none(getattr(options, "fq2", None))

    total_variants, variants_by_contig = _count_variants(vcf_path)
    total_reads = _count_reads(bam_path, fq1_path, fq2_path, bool(getattr(options, "paired_ended", False)))

    fastq_outputs = [p for p in (fq1_path, fq2_path) if p is not None]

    summary = {
        "schema_version": SCHEMA_VERSION,
        "neat_version": NEAT_VERSION,
        "run": {
            "started_at": started_at,
            "completed_at": completed_at,
            "duration_seconds": round(completed_ts - analysis_start, 3),
            "config_file": _abs_or_none(config_path),
            "output_dir": str(Path(output_dir).resolve()),
            "output_prefix": file_prefix,
        },
        "config": {
            "reference":       _abs_or_none(getattr(options, "reference", None)),
            "coverage":        getattr(options, "coverage", None),
            "read_len":        getattr(options, "read_len", None),
            "paired_ended":    getattr(options, "paired_ended", None),
            "fragment_mean":   getattr(options, "fragment_mean", None),
            "fragment_st_dev": getattr(options, "fragment_st_dev", None),
            "ploidy":          getattr(options, "ploidy", None),
            "rng_seed":        getattr(options, "rng_seed", None),
            "threads":         getattr(options, "threads", None),
            "mutation_rate":   getattr(options, "mutation_rate", None),
            "mutation_bed":    _abs_or_none(getattr(options, "mutation_bed", None)),
            "target_bed":      _abs_or_none(getattr(options, "target_bed", None)),
            "discard_bed":     _abs_or_none(getattr(options, "discard_bed", None)),
            "include_vcf":     _abs_or_none(getattr(options, "include_vcf", None)),
            "mutation_model":  _abs_or_none(getattr(options, "mutation_model", None)),
            "gc_model":        _abs_or_none(getattr(options, "gc_model", None)),
            "error_model":     _abs_or_none(getattr(options, "error_model", None)),
            "fragment_model":  _abs_or_none(getattr(options, "fragment_model", None)),
        },
        "outputs": {
            "fastq": fastq_outputs if fastq_outputs else None,
            "bam":   bam_path,
            "vcf":   vcf_path,
        },
        "delivered": {
            "total_reads":        total_reads,
            "total_variants":     total_variants,
            "variants_by_contig": variants_by_contig,
            "contigs_simulated":  list(contigs_simulated),
            "reference_contigs":  list(reference_contigs) if reference_contigs is not None else list(contigs_simulated),
        },
    }

    out_path = Path(output_dir) / "simulation_summary.json"
    tmp_path = out_path.with_suffix(".json.tmp")
    with open(tmp_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    os.replace(tmp_path, out_path)
    _LOG.info(f"Wrote {out_path}")
    return out_path


def _iso_utc(epoch_seconds: float) -> str:
    return datetime.fromtimestamp(epoch_seconds, tz=timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _abs_or_none(value):
    if value is None:
        return None
    return str(Path(value).resolve())


def _count_variants(vcf_path):
    if vcf_path is None or not Path(vcf_path).is_file():
        return None, None
    total = 0
    by_contig: dict[str, int] = {}
    try:
        with pysam.VariantFile(vcf_path) as vf:
            for rec in vf:
                total += 1
                by_contig[rec.chrom] = by_contig.get(rec.chrom, 0) + 1
    except Exception as exc:
        _LOG.warning(f"Could not count variants in {vcf_path}: {exc}")
        return None, None
    return total, by_contig


def _count_reads(bam_path, fq1_path, fq2_path, paired_ended):
    if bam_path is not None and Path(bam_path).is_file():
        try:
            with pysam.AlignmentFile(bam_path, "rb") as bf:
                return bf.count(until_eof=True)
        except Exception as exc:
            _LOG.warning(f"Could not count reads in {bam_path}: {exc}")
    if fq1_path is not None and Path(fq1_path).is_file():
        try:
            n = sum(1 for _ in pysam.FastxFile(fq1_path))
            return n * 2 if paired_ended else n
        except Exception as exc:
            _LOG.warning(f"Could not count reads in {fq1_path}: {exc}")
    return None
