"""
Run the simulator across a configuration matrix and assert the output describes itself.

The axes here are the ones that have actually surfaced defects, not every option NEAT has:

  read_len              only the 250 bp configuration exposed #335; 75 and 150 were clean
  fragment_mean vs      short inserts make the mates cover the same window, which is what
    read_len            broke BAM sort order, and remove the headroom a deletion needs
  adapters              soft-clip placement is strand-dependent and easy to get backwards
  paired vs single      single-ended never produces a reverse read, so it exercises none of
                        the reverse-strand coordinate handling
  include_vcf           input variants reach the golden VCF through a different path than
                        simulated ones, and MNPs become UnknownVariant
  reference composition repeats and skewed GC change what the aligner in the CIGAR fallback
                        has to work with
"""

from __future__ import annotations

import pytest

from .conftest import requires_integration
from . import invariants

pytestmark = requires_integration


# (id, reference key, config overrides)
MATRIX = [
    ("plain_150",        "uniform",      {"read_len": 150, "coverage": 20, "paired_ended": True,
                                          "fragment_mean": 350, "fragment_st_dev": 50,
                                          "mutation_rate": 0.001}),
    ("plain_250",        "uniform",      {"read_len": 250, "coverage": 20, "paired_ended": True,
                                          "fragment_mean": 350, "fragment_st_dev": 50,
                                          "mutation_rate": 0.001}),
    ("plain_75",         "uniform",      {"read_len": 75, "coverage": 20, "paired_ended": True,
                                          "fragment_mean": 300, "fragment_st_dev": 50,
                                          "mutation_rate": 0.001}),
    ("high_mutation",    "uniform",      {"read_len": 150, "coverage": 15, "paired_ended": True,
                                          "fragment_mean": 350, "fragment_st_dev": 50,
                                          "mutation_rate": 0.01}),
    ("single_ended",     "uniform",      {"read_len": 150, "coverage": 15, "paired_ended": False,
                                          "fragment_mean": 350, "fragment_st_dev": 50,
                                          "mutation_rate": 0.001}),
    ("gc_rich",          "gc_rich",      {"read_len": 150, "coverage": 15, "paired_ended": True,
                                          "fragment_mean": 350, "fragment_st_dev": 50,
                                          "mutation_rate": 0.001}),
    ("repetitive",       "repetitive",   {"read_len": 150, "coverage": 15, "paired_ended": True,
                                          "fragment_mean": 350, "fragment_st_dev": 50,
                                          "mutation_rate": 0.001}),
    ("multi_contig",     "multi_contig", {"read_len": 150, "coverage": 15, "paired_ended": True,
                                          "fragment_mean": 350, "fragment_st_dev": 50,
                                          "mutation_rate": 0.001}),
    ("adapters_truseq",  "uniform",      {"read_len": 150, "coverage": 15, "paired_ended": True,
                                          "fragment_mean": 130, "fragment_st_dev": 40,
                                          "adapters": True, "adapter_preset": "truseq",
                                          "mutation_rate": 0.001}),
    ("adapters_nextera", "uniform",      {"read_len": 150, "coverage": 15, "paired_ended": True,
                                          "fragment_mean": 130, "fragment_st_dev": 40,
                                          "adapters": True, "adapter_preset": "nextera",
                                          "mutation_rate": 0.001}),
    ("adapters_single",  "uniform",      {"read_len": 150, "coverage": 15, "paired_ended": False,
                                          "fragment_mean": 130, "fragment_st_dev": 40,
                                          "adapters": True, "adapter_preset": "truseq",
                                          "mutation_rate": 0.001}),
    ("adapters_saturated", "uniform",    {"read_len": 150, "coverage": 15, "paired_ended": True,
                                          "fragment_mean": 60, "fragment_st_dev": 15,
                                          "adapters": True, "adapter_preset": "truseq",
                                          "mutation_rate": 0.0}),
    ("keep_short_only",  "uniform",      {"read_len": 150, "coverage": 15, "paired_ended": True,
                                          "fragment_mean": 130, "fragment_st_dev": 40,
                                          "keep_short_fragments": True, "mutation_rate": 0.001}),
]

ADAPTER_CONFIGS = {"adapters_truseq", "adapters_nextera", "adapters_single", "adapters_saturated"}

# Configurations with a known, open defect. Recorded here so the suite is green when nothing
# *new* is wrong — a nightly that is permanently red gets ignored — while still failing if the
# defect spreads to a configuration not listed. Delete an entry when its issue closes; the test
# then simply passes.
KNOWN_TRAILING_DELETION = {
    "keep_short_only": (
        "302eaa1 deliberately keeps a deletion that reaches a short insert's last base, since "
        "there is no headroom past the fragment to absorb it. SAM has no way to express that: a "
        "CIGAR ending in D describes reference past the read's last base, overstating the span, "
        "and Picard's ValidateSamFile rejects it. Representable only by shortening the alignment."
    ),
}

KNOWN_MISPLACED = {
    "plain_250": (
        "#335: errors are applied in the order the sampler accepted them rather than in position "
        "order, so an error behind the high-water mark lengthens the read without consuming a "
        "base. Fixing it needs the quality accounting in #337 untangled first."
    ),
}


@pytest.fixture(scope="module")
def _runs():
    """Simulations are the expensive part, so each configuration is run once and its output
    shared by every check below."""
    return {}


def _run_for(case_id, reference_key, config, references, run_simulation, _runs):
    if case_id not in _runs:
        _runs[case_id] = run_simulation(references[reference_key], config, prefix=case_id)
    return _runs[case_id]


@pytest.mark.parametrize("case_id,reference_key,config", MATRIX, ids=[c[0] for c in MATRIX])
def test_reads_are_placed_where_their_alignment_says(
    case_id, reference_key, config, references, run_simulation, _runs
):
    run = _run_for(case_id, reference_key, config, references, run_simulation, _runs)
    failures = invariants.reads_are_placed_where_they_say(run.bam, run.reference)
    if failures and case_id in KNOWN_MISPLACED:
        pytest.xfail(f"{KNOWN_MISPLACED[case_id]} ({len(failures)} records)")
    assert not failures, (
        f"[{case_id}] {len(failures)} misplaced record(s):\n  "
        + "\n  ".join(failures[:10])
    )


@pytest.mark.parametrize("case_id,reference_key,config", MATRIX, ids=[c[0] for c in MATRIX])
def test_cigars_are_consistent(case_id, reference_key, config, references, run_simulation, _runs):
    run = _run_for(case_id, reference_key, config, references, run_simulation, _runs)
    failures = invariants.cigars_account_for_every_base(run.bam)
    well_formed = invariants.cigars_are_well_formed(run.bam)
    if well_formed and case_id in KNOWN_TRAILING_DELETION:
        if not failures:
            pytest.xfail(f"{KNOWN_TRAILING_DELETION[case_id]} ({len(well_formed)} records)")
    else:
        failures += well_formed
    assert not failures, f"[{case_id}] {len(failures)} bad cigar(s):\n  " + "\n  ".join(failures[:10])


@pytest.mark.parametrize("case_id,reference_key,config", MATRIX, ids=[c[0] for c in MATRIX])
def test_golden_bam_is_sorted_and_indexable(
    case_id, reference_key, config, references, run_simulation, _runs
):
    run = _run_for(case_id, reference_key, config, references, run_simulation, _runs)
    failures = invariants.bam_is_sorted_and_indexable(run.bam)
    assert not failures, f"[{case_id}]:\n  " + "\n  ".join(failures[:10])


@pytest.mark.parametrize("case_id,reference_key,config", MATRIX, ids=[c[0] for c in MATRIX])
def test_fastq_and_bam_agree(case_id, reference_key, config, references, run_simulation, _runs):
    run = _run_for(case_id, reference_key, config, references, run_simulation, _runs)
    failures = invariants.fastq_records_match_the_bam(run.bam, run.fastqs)
    assert not failures, f"[{case_id}]:\n  " + "\n  ".join(failures[:10])


@pytest.mark.parametrize(
    "case_id,reference_key,config",
    [c for c in MATRIX if c[0] in ADAPTER_CONFIGS],
    ids=[c[0] for c in MATRIX if c[0] in ADAPTER_CONFIGS],
)
def test_adapter_clips_are_on_the_right_side(
    case_id, reference_key, config, references, run_simulation, _runs
):
    run = _run_for(case_id, reference_key, config, references, run_simulation, _runs)
    failures = invariants.adapter_clips_are_on_the_right_side(run.bam)
    assert not failures, f"[{case_id}]:\n  " + "\n  ".join(failures[:10])
