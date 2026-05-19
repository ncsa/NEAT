"""
NEAT-aware false-negative attribution for `neat compare-vcfs`.

Each FN from hap.py is tagged with one or more reasons drawn from the
simulator's own configuration:

  - `outside_simulated_contigs` — the FN's chromosome wasn't in the NEAT run.
  - `outside_mutation_bed`      — `mutation_bed` was set and the FN position
                                  falls outside its regions.
  - `outside_target_bed`        — `target_bed` was set and the FN position
                                  falls outside its regions.
  - `unknown`                   — none of the above; NEAT can't explain it.

If the FN's contig wasn't simulated at all, that single tag is reported on its
own (the BED checks would be meaningless).
"""
import logging
from pathlib import Path

__all__ = [
    "attribute_fn",
    "attribute_fns",
    "load_bed_intervals",
    "position_in_intervals",
    "REASON_OUTSIDE_CONTIGS",
    "REASON_OUTSIDE_MUTATION_BED",
    "REASON_OUTSIDE_TARGET_BED",
    "REASON_UNKNOWN",
]

_LOG = logging.getLogger(__name__)

REASON_OUTSIDE_CONTIGS = "outside_simulated_contigs"
REASON_OUTSIDE_MUTATION_BED = "outside_mutation_bed"
REASON_OUTSIDE_TARGET_BED = "outside_target_bed"
REASON_UNKNOWN = "unknown"


def load_bed_intervals(bed_path: Path | str | None) -> dict[str, list[tuple[int, int]]] | None:
    """
    Parse a BED file into a per-contig list of sorted 0-based half-open
    intervals. Returns None if `bed_path` is None.

    Comments (`#`), `track`, and `browser` header lines are skipped. Rows with
    non-integer start/end are dropped with a debug log entry.
    """
    if bed_path is None:
        return None
    intervals: dict[str, list[tuple[int, int]]] = {}
    with open(bed_path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith(("track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                _LOG.debug(f"{bed_path}:{lineno}: skipping short line")
                continue
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                _LOG.debug(f"{bed_path}:{lineno}: skipping non-integer interval")
                continue
            intervals.setdefault(parts[0], []).append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def position_in_intervals(pos: int, intervals: list[tuple[int, int]]) -> bool:
    """
    Membership test for a 1-based VCF position against 0-based half-open BED
    intervals sorted by start. Returns True if any interval contains the
    position. Correct for overlapping intervals.
    """
    pos_0 = pos - 1
    for start, end in intervals:
        if start > pos_0:
            return False  # sorted: no further interval can match
        if pos_0 < end:
            return True
    return False


def attribute_fn(
    chrom: str,
    pos: int,
    contigs_simulated: set[str] | frozenset[str],
    mutation_intervals: dict[str, list[tuple[int, int]]] | None,
    target_intervals: dict[str, list[tuple[int, int]]] | None,
) -> list[str]:
    """
    Return the list of NEAT-specific reasons that explain a single FN.

    If the contig wasn't simulated, that root cause is reported alone — BED
    checks are skipped because they presuppose simulation.
    """
    if chrom not in contigs_simulated:
        return [REASON_OUTSIDE_CONTIGS]

    reasons: list[str] = []
    if mutation_intervals is not None:
        chrom_intervals = mutation_intervals.get(chrom, [])
        if not position_in_intervals(pos, chrom_intervals):
            reasons.append(REASON_OUTSIDE_MUTATION_BED)
    if target_intervals is not None:
        chrom_intervals = target_intervals.get(chrom, [])
        if not position_in_intervals(pos, chrom_intervals):
            reasons.append(REASON_OUTSIDE_TARGET_BED)
    if not reasons:
        reasons.append(REASON_UNKNOWN)
    return reasons


def attribute_fns(fn_records, summary: dict) -> list[tuple]:
    """
    Tag every FN against the run's simulation_summary.

    :param fn_records: iterable of pysam.VariantRecord (FN bucket from hap.py).
    :param summary: parsed simulation_summary.json.
    :return: list of (record, reasons) tuples; `reasons` is a list[str].
    """
    contigs = frozenset(summary["delivered"].get("contigs_simulated", []))
    cfg = summary.get("config", {})
    mutation_intervals = load_bed_intervals(cfg.get("mutation_bed"))
    target_intervals = load_bed_intervals(cfg.get("target_bed"))

    return [
        (rec, attribute_fn(rec.chrom, rec.pos, contigs, mutation_intervals, target_intervals))
        for rec in fn_records
    ]
