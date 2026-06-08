"""
hap.py subprocess invocation + output-VCF parsing for `neat compare-vcfs`.

hap.py emits one VCF per run with per-sample FORMAT annotations describing each
record's classification:
  - sample 0 (TRUTH): BD == 'TP' (matched), 'FN' (missed), or '.' (no call)
  - sample 1 (QUERY): BD == 'TP' (matched), 'FP' (spurious), or '.' (no call)

This module wraps the subprocess and classifies each record into one of TP/FN/FP
based on those FORMAT fields. Multi-allelic and gVCF handling are delegated to
hap.py upstream; we read what it emits.
"""
import logging
import subprocess
from pathlib import Path
from typing import Iterable

import pysam

__all__ = [
    "HappyExecutionError",
    "HappyParseError",
    "run_happy",
    "parse_happy_output",
]

_LOG = logging.getLogger(__name__)


class HappyExecutionError(RuntimeError):
    """Raised when the hap.py subprocess fails."""


class HappyParseError(RuntimeError):
    """Raised when the hap.py output VCF is missing expected structure."""


def run_happy(
    happy_bin: Path,
    golden_vcf: Path,
    called_vcf: Path,
    output_prefix: Path,
    reference: Path | None = None,
    target_bed: Path | None = None,
    extra_args: Iterable[str] = (),
) -> Path:
    """
    Invoke hap.py and return the path to its bgzipped output VCF.

    :param happy_bin: Absolute path to the hap.py executable.
    :param golden_vcf: Truth VCF (NEAT golden).
    :param called_vcf: Query VCF (downstream variant caller).
    :param output_prefix: Path prefix; hap.py appends `.vcf.gz` and writes
        sibling artifacts (summary.csv, extended.csv, runinfo.json).
    :param reference: Optional reference FASTA forwarded as `-r`.
    :param target_bed: Optional restricting BED forwarded as `-T`.
    :param extra_args: Optional positional pass-through for future tunables.
    :return: Path to `<output_prefix>.vcf.gz`.
    :raises HappyExecutionError: if hap.py exits non-zero or its output VCF is absent.
    """
    cmd = [
        str(happy_bin), str(golden_vcf), str(called_vcf),
        "-o", str(output_prefix),
    ]
    if reference is not None:
        cmd += ["-r", str(reference)]
    if target_bed is not None:
        cmd += ["-T", str(target_bed)]
    cmd += list(extra_args)

    _LOG.info(f"Running hap.py: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        _LOG.error(f"hap.py stderr:\n{result.stderr}")
        raise HappyExecutionError(
            f"hap.py exited {result.returncode}. See log for stderr. "
            f"Command: {' '.join(cmd)}"
        )
    if result.stderr:
        _LOG.debug(f"hap.py stderr:\n{result.stderr}")

    output_vcf = Path(str(output_prefix) + ".vcf.gz")
    if not output_vcf.is_file():
        raise HappyExecutionError(
            f"hap.py reported success but its output VCF is missing: {output_vcf}"
        )
    return output_vcf


def parse_happy_output(vcf_path: Path) -> dict[str, list]:
    """
    Group hap.py output records into TP/FN/FP buckets.

    :param vcf_path: hap.py's bgzipped output VCF.
    :return: dict with keys 'TP', 'FN', 'FP'; values are lists of pysam.VariantRecord.
    :raises HappyParseError: if the VCF lacks the expected TRUTH/QUERY samples
        or the BD FORMAT field.
    """
    buckets: dict[str, list] = {"TP": [], "FN": [], "FP": []}
    with pysam.VariantFile(str(vcf_path)) as vf:
        sample_names = list(vf.header.samples)
        if len(sample_names) < 2:
            raise HappyParseError(
                f"{vcf_path} has {len(sample_names)} sample(s); hap.py output must have "
                f"TRUTH and QUERY samples."
            )
        if "BD" not in vf.header.formats:
            raise HappyParseError(
                f"{vcf_path} is missing the BD FORMAT field — not a hap.py output VCF?"
            )
        truth_name, query_name = sample_names[0], sample_names[1]

        for rec in vf:
            truth_bd = rec.samples[truth_name].get("BD") or "."
            query_bd = rec.samples[query_name].get("BD") or "."
            classification = _classify(truth_bd, query_bd)
            if classification is not None:
                buckets[classification].append(rec)
    return buckets


def _classify(truth_bd: str, query_bd: str) -> str | None:
    """
    Apply the hap.py decision rules.

    A record is:
      - TP if either sample reports TP (matched on at least one side)
      - FN if truth reports FN and query did not match
      - FP if query reports FP and truth did not match
      - Otherwise None (no-call / nocomp / hap.py-internal types)
    """
    if truth_bd == "TP" or query_bd == "TP":
        return "TP"
    if truth_bd == "FN":
        return "FN"
    if query_bd == "FP":
        return "FP"
    return None
