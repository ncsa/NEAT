"""
Properties every NEAT run's output must satisfy, checked over every record it produced.

The unit suite pins behaviour that is already understood: a known bug, a known branch. What it
cannot do is notice that 85 of 16,000 reads are one base off, because the run exits 0, writes a
complete output set, and every record is individually well formed. Every defect these checks were
written for (#324, #325, #326, #335) shipped past a green suite and was found by generating reads
and asking whether the output describes itself consistently.

Each check returns a list of human-readable failure strings rather than asserting, so one run can
report everything that is wrong with it instead of stopping at the first record.
"""

from __future__ import annotations

import gzip
import re
import subprocess
from pathlib import Path

import pysam
from Bio import SeqIO

# A record is judged misplaced only if it matches the reference substantially better somewhere
# else. Sequencing error alone can push the mismatch rate high — on a 250 bp read the default
# model's tail runs near 18% — so an absolute threshold cannot separate "noisy" from "misplaced",
# but a large relative improvement at a nearby offset can.
_MISPLACED_FLOOR = 0.30
_MISPLACED_RATIO = 3.0
_SHIFTS = (-3, -2, -1, 1, 2, 3)


def _load_reference(reference_path: Path) -> dict[str, str]:
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(reference_path), "fasta")}


def _mismatch_rate(record, reference: str, seq: str, start: int) -> float:
    """Mismatch rate over the record's M columns, walking its CIGAR from `start`."""
    ref_pos, query_pos, bad, total = start, 0, 0, 0
    for op, length in record.cigartuples:
        if op == 0:                                  # M
            window = reference[ref_pos:ref_pos + length]
            chunk = seq[query_pos:query_pos + length]
            bad += sum(1 for a, b in zip(chunk, window) if a != b)
            total += min(len(window), len(chunk))
            ref_pos += length
            query_pos += length
        elif op == 1:                                # I — query only
            query_pos += length
        elif op == 2:                                # D — reference only
            ref_pos += length
        elif op == 4:                                # S — soft clip, query only
            query_pos += length
    return bad / total if total else 1.0


def reads_are_placed_where_they_say(bam_path: Path, reference_path: Path) -> list[str]:
    """
    Every record's POS and CIGAR must describe where its sequence actually came from.

    This is the check that found #326 (gapped reverse reads walked backwards through the aligned
    columns) and #335 (a duplicated error position lengthening the read). Both produced records
    that claimed a clean alignment at the wrong coordinate.
    """
    reference = _load_reference(reference_path)
    failures = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for record in bam:
            if record.is_unmapped or not record.cigartuples or not record.query_sequence:
                continue
            contig = reference.get(bam.get_reference_name(record.reference_id))
            if contig is None:
                failures.append(f"{record.query_name}: contig not in reference")
                continue
            seq = record.query_sequence.upper()
            here = _mismatch_rate(record, contig, seq, record.reference_start)
            if here <= _MISPLACED_FLOOR:
                continue
            for shift in _SHIFTS:
                start = record.reference_start + shift
                if start < 0:
                    continue
                there = _mismatch_rate(record, contig, seq, start)
                if there * _MISPLACED_RATIO < here and there < 0.10:
                    failures.append(
                        f"{record.query_name}: POS {record.reference_start} {record.cigarstring} "
                        f"mismatch {here:.0%}, but {there:.0%} at POS{shift:+d}"
                    )
                    break
    return failures


def cigars_account_for_every_base(bam_path: Path) -> list[str]:
    """M + I + S must cover the whole read: those are the query-consuming ops.

    A CIGAR that comes up short means bases went somewhere unrecorded — #326's insertion handling
    dropped everything past a fixed-size op list this way.
    """
    failures = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for record in bam:
            if not record.cigartuples or not record.query_sequence:
                continue
            covered = sum(n for op, n in record.cigartuples if op in (0, 1, 4))
            if covered != len(record.query_sequence):
                failures.append(
                    f"{record.query_name}: cigar {record.cigarstring} covers {covered} bases "
                    f"but SEQ is {len(record.query_sequence)}"
                )
    return failures


def cigars_have_no_leading_deletion(bam_path: Path) -> list[str]:
    """A CIGAR cannot open on a deletion: it would consume reference before the read's first
    base, which is what POS already fixes. Nothing in NEAT is known to produce one, so this is a
    hard failure wherever it appears."""
    failures = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for record in bam:
            ops = record.cigartuples
            if ops and ops[0][0] == 2:
                failures.append(f"{record.query_name}: cigar {record.cigarstring} opens on D")
    return failures


def cigars_have_no_trailing_deletion(bam_path: Path) -> list[str]:
    """A CIGAR ending in a deletion describes reference past the read's last base, overstating
    the span, and Picard's ValidateSamFile rejects it.

    Kept separate from the leading case because this one is a known, deliberate tradeoff rather
    than an unexpected defect — see the known-failure note in test_output_invariants.py.
    """
    failures = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for record in bam:
            ops = record.cigartuples
            if ops and ops[-1][0] == 2:
                failures.append(f"{record.query_name}: cigar {record.cigarstring} ends on D")
    return failures


def bam_is_sorted_and_indexable(bam_path: Path) -> list[str]:
    """
    Coordinate order, then an actual index build.

    Sortedness is not cosmetic: the runner indexes the golden BAM, so an out-of-order record
    aborts the whole run after the FASTQs are already on disk. Checking the order rather than only
    the index gives the offending pair by name.
    """
    failures = []
    previous = None
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for record in bam:
            if record.is_unmapped or not record.cigartuples:
                continue
            # Coordinate order is per contig, and position resets at each contig boundary, so
            # compare within a contig and require the contig index itself to be non-decreasing.
            current = (record.reference_id, record.reference_start)
            if previous is not None:
                if current[0] < previous[1][0]:
                    failures.append(
                        f"contig out of order: {record.query_name} on ref {current[0]} "
                        f"follows {previous[0]} on ref {previous[1][0]}"
                    )
                elif current[0] == previous[1][0] and current[1] < previous[1][1]:
                    failures.append(
                        f"out of order: {record.query_name} at {current[1]} "
                        f"follows {previous[0]} at {previous[1][1]}"
                    )
            previous = (record.query_name, current)
    try:
        pysam.index(str(bam_path))
    except Exception as exc:                                     # pragma: no cover - env dependent
        failures.append(f"samtools index failed: {exc}")
    return failures


def fastq_records_match_the_bam(bam_path: Path, fastq_paths: list[Path]) -> list[str]:
    """Same number of reads in the FASTQs as in the golden BAM, and one quality character per
    base in every FASTQ record."""
    failures = []
    fastq_reads = 0
    for path in fastq_paths:
        with gzip.open(path, "rt") as handle:
            for index, line in enumerate(handle):
                if index % 4 == 1:
                    fastq_reads += 1
                    bases = len(line.strip())
                elif index % 4 == 3:
                    scores = len(line.strip())
                    if scores != bases:
                        failures.append(f"{path.name}: {scores} quality scores for {bases} bases")
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        bam_reads = sum(1 for _ in bam)
    if fastq_reads != bam_reads:
        failures.append(f"{fastq_reads} FASTQ reads but {bam_reads} BAM records")
    return failures


def vcf_is_valid(vcf_path: Path) -> list[str]:
    """
    bcftools has to accept the golden VCF.

    NEAT itself runs `bcftools sort` over this file, so a record referencing an undeclared INFO
    key or FILTER value does not merely produce an odd file — it fails the run at the very end
    (#325). Reading it back is the cheapest way to assert the header and records agree.
    """
    result = subprocess.run(
        ["bcftools", "view", "-h", str(vcf_path)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        return [f"bcftools rejected {vcf_path.name}: {result.stderr.strip()[:200]}"]

    header = result.stdout
    declared = set(re.findall(r"##(?:INFO|FILTER|FORMAT)=<ID=([^,>]+)", header))
    failures = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for record in vcf:
            for key in record.info:
                if key not in declared:
                    failures.append(f"{record.chrom}:{record.pos} uses undeclared INFO key {key}")
            for value in (record.filter.keys() or []):
                if value not in declared and value != "PASS":
                    failures.append(f"{record.chrom}:{record.pos} uses undeclared FILTER {value}")
    return failures[:20]


def adapter_clips_are_on_the_right_side(bam_path: Path) -> list[str]:
    """
    SEQ is stored reference-forward, and the adapter sits at the read's 3' end — so the soft clip
    belongs at the start of the CIGAR for a reverse read and at the end for a forward one. Getting
    it backwards still passes `samtools quickcheck`; it just misplaces every short-insert read.
    """
    failures = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for record in bam:
            ops = record.cigartuples
            if not ops:
                continue
            if record.is_reverse and ops[-1][0] == 4:
                failures.append(f"{record.query_name}: reverse read soft-clipped at its end")
            if not record.is_reverse and ops[0][0] == 4:
                failures.append(f"{record.query_name}: forward read soft-clipped at its start")
    return failures[:20]
