"""
Chromosome-name normalization utilities for cross-checking reference/VCF/BED
conventions.

Genomics tooling ships with multiple incompatible naming conventions:
  - UCSC:    `chr1`, `chr2`, ..., `chrX`, `chrM`
  - Ensembl: `1`,    `2`,    ..., `X`,    `MT`
  - NCBI accession: `NC_000001.11`, ...
  - Custom:  whatever the lab named the scaffolds

A user can run NEAT with a reference using one convention but supply a
`mutation_bed` / `target_bed` using another. NEAT does not silently normalize;
instead, callers detect a mismatch up front and either warn (see
`compare_vcfs`) or apply an explicit user-supplied alias map.

This module is the single source of truth for chrom-name heuristics so the
simulator side and the compare-vcfs side stay consistent.
"""
import logging
from pathlib import Path

__all__ = [
    "apply_aliases",
    "find_aliases",
    "load_chrom_aliases",
    "prefix_flip_candidates",
]

_LOG = logging.getLogger(__name__)

# Common mitochondrial name variants. A reference using any of these
# represents the same physical chromosome.
_MITO_VARIANTS = frozenset({"M", "MT", "chrM", "chrMT"})


def prefix_flip_candidates(name: str) -> set[str]:
    """
    Return the set of alternative names that probably refer to the same
    contig as `name`. Excludes `name` itself.

    Two heuristics:
      1. `chr` prefix add/strip — `chr1` ↔ `1`, `chrX` ↔ `X`.
      2. Mitochondrial variants — names in {M, MT, chrM, chrMT} are all
         treated as aliases of each other.

    No other heuristic — accession-style names and custom scaffold names
    require an explicit user-supplied alias map.
    """
    candidates: set[str] = set()
    if name.startswith("chr"):
        candidates.add(name[3:])
    else:
        candidates.add(f"chr{name}")
    if name in _MITO_VARIANTS:
        candidates |= _MITO_VARIANTS
    candidates.discard(name)
    return candidates


def find_aliases(source: set[str] | frozenset[str], target: set[str] | frozenset[str]) -> dict[str, str]:
    """
    For each name in `source` that doesn't appear natively in `target`, propose
    a canonical form from `target` if one of the heuristics matches. Returns a
    `{source_name: target_name}` dict; names with a native match or no inferrable
    candidate are omitted.
    """
    mapping: dict[str, str] = {}
    for name in source:
        if name in target:
            continue
        for candidate in prefix_flip_candidates(name):
            if candidate in target:
                mapping[name] = candidate
                break
    return mapping


def load_chrom_aliases(path: Path | str | None) -> dict[str, str]:
    """
    Parse a two-column alias TSV: `source_name<TAB>canonical_name`. Lines
    starting with `#` and empty lines are skipped. Returns `{}` when `path`
    is None.

    Whitespace tolerated as a fallback separator for human-typed files.
    """
    if path is None:
        return {}
    aliases: dict[str, str] = {}
    with open(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                _LOG.warning(f"{path}:{lineno}: skipping malformed alias line: {raw!r}")
                continue
            aliases[parts[0]] = parts[1]
    return aliases


def apply_aliases(name: str, aliases: dict[str, str]) -> str:
    """Return the canonical form of `name`, or `name` unchanged if no alias maps it."""
    return aliases.get(name, name)
