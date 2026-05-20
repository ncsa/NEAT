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

from ..common.chrom_names import find_aliases

__all__ = [
    "attribute_fn",
    "attribute_fns",
    "detect_chrom_naming_mismatches",
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


def load_bed_intervals(
    bed_path: Path | str | None,
    aliases: dict[str, str] | None = None,
) -> dict[str, list[tuple[int, int]]] | None:
    """
    Parse a BED file into a per-contig list of sorted 0-based half-open
    intervals. Returns None if `bed_path` is None.

    Comments (`#`), `track`, and `browser` header lines are skipped. Rows with
    non-integer start/end are dropped with a debug log entry.

    If `aliases` is provided, each BED row's chrom field is remapped through
    the dict (`aliases.get(chrom, chrom)`) before being stored, so downstream
    lookups can use reference-canonical names.
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
            chrom = parts[0]
            if aliases:
                chrom = aliases.get(chrom, chrom)
            intervals.setdefault(chrom, []).append((start, end))
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


def attribute_fns(
    fn_records,
    summary: dict,
    aliases: dict[str, str] | None = None,
    skip_beds: set[str] | None = None,
) -> list[tuple]:
    """
    Tag every FN against the run's simulation_summary.

    :param fn_records: iterable of pysam.VariantRecord (FN bucket from hap.py).
    :param summary: parsed simulation_summary.json.
    :param aliases: optional user-supplied {bed_name: canonical_name} map applied
        to BED chrom names at load time.
    :param skip_beds: optional set of BED labels (e.g., {"mutation_bed"}) to
        treat as if not configured. Used by the runner to skip BEDs whose chrom
        names are entirely mismatched against the reference — attributing FNs
        against an unusable BED would produce misleading `outside_*` counts.
    :return: list of (record, reasons) tuples; `reasons` is a list[str].
    """
    skip_beds = skip_beds or set()
    contigs = frozenset(summary["delivered"].get("contigs_simulated", []))
    cfg = summary.get("config", {})
    mutation_intervals = (
        None if "mutation_bed" in skip_beds
        else load_bed_intervals(cfg.get("mutation_bed"), aliases=aliases)
    )
    target_intervals = (
        None if "target_bed" in skip_beds
        else load_bed_intervals(cfg.get("target_bed"), aliases=aliases)
    )

    return [
        (rec, attribute_fn(rec.chrom, rec.pos, contigs, mutation_intervals, target_intervals))
        for rec in fn_records
    ]


def detect_chrom_naming_mismatches(
    summary: dict,
    aliases: dict[str, str] | None = None,
) -> list[dict]:
    """
    Inspect each configured BED's chrom set against the reference's; return a
    list of warning records (one per mismatched BED) or empty if all BEDs
    overlap (after any user-supplied aliases are applied).

    The "reference" set comes from `summary['delivered']['reference_contigs']`
    (the full FASTA contig list), with a fallback to `contigs_simulated` for
    backward compatibility with summaries written before that field existed.
    """
    aliases = aliases or {}
    reference_chroms = frozenset(
        summary["delivered"].get("reference_contigs")
        or summary["delivered"].get("contigs_simulated", [])
    )
    if not reference_chroms:
        return []

    warnings: list[dict] = []
    for bed_label in ("mutation_bed", "target_bed"):
        bed_path = summary.get("config", {}).get(bed_label)
        if not bed_path:
            continue
        try:
            intervals = load_bed_intervals(bed_path, aliases=None)
        except FileNotFoundError:
            continue  # surface this separately; not a naming issue
        if not intervals:
            continue

        raw_chroms = frozenset(intervals.keys())
        mapped_chroms = frozenset(aliases.get(c, c) for c in raw_chroms)
        if mapped_chroms & reference_chroms:
            continue  # at least one chrom matches; not a mismatch worth flagging

        suggested = find_aliases(raw_chroms, reference_chroms)
        message = (
            f"{bed_label} chrom names don't overlap with reference contigs"
            + (
                f"; suggested --chrom-aliases mappings: {suggested}"
                if suggested
                else " and no naming convention could be inferred. Attribution against "
                "this BED will report every FN as 'outside'."
            )
        )
        warnings.append({
            "type": "chrom_naming_mismatch",
            "bed": bed_label,
            "bed_path": str(bed_path),
            "bed_chroms_sample": sorted(raw_chroms)[:5],
            "reference_chroms_sample": sorted(reference_chroms)[:5],
            "suggested_aliases": suggested,
            "message": message,
        })
    return warnings
