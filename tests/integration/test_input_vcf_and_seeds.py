"""
Two things the configuration matrix cannot cover on its own.

The first is the input-VCF path: variants read from `include_vcf` reach the golden VCF through
different code than simulated ones, and their metadata has to survive the trip. #325 broke that
for every input variant at once and no test noticed, because nothing read the output back.

The second is seed dependence. An invariant that holds for one seed can be holding by
cancellation rather than by construction — #337 is exactly that, a quality-length accounting
that balances only in the order one particular draw produces. Checking across seeds is what
tells the two apart.
"""

from __future__ import annotations

import pytest

from .conftest import requires_integration
from . import invariants

pytestmark = requires_integration


@pytest.fixture(scope="module")
def input_vcf_run(references, mixed_variant_vcf, run_simulation):
    return run_simulation(
        references["uniform"],
        {
            "read_len": 150, "coverage": 20, "paired_ended": True,
            "fragment_mean": 350, "fragment_st_dev": 50, "mutation_rate": 0.001,
            "include_vcf": str(mixed_variant_vcf), "produce_vcf": True,
        },
        prefix="input_vcf",
    )


def test_input_variants_survive_into_a_valid_golden_vcf(input_vcf_run):
    """Every INFO key and FILTER value a record uses must be declared in the header NEAT wrote.

    NEAT runs `bcftools sort` over this file itself, so an undeclared key does not just make an
    odd VCF — it fails the run after the FASTQs are on disk (#325).
    """
    failures = invariants.vcf_is_valid(input_vcf_run.vcf)
    assert not failures, "golden VCF is not self-consistent:\n  " + "\n  ".join(failures)


def test_every_input_variant_reaches_the_output(input_vcf_run, mixed_variant_vcf):
    """The no-call record is the interesting one: it has to be kept with a generated genotype
    rather than dropped (#324) or crashed on."""
    import pysam

    wanted = {"rs_snv", "rs_nocall", "rs_del", "rs_ins", "rs_mnp", "rs_multi"}
    with pysam.VariantFile(str(input_vcf_run.vcf)) as vcf:
        found = {record.id for record in vcf if record.id}
    missing = wanted - found
    assert not missing, f"input variants missing from the golden VCF: {sorted(missing)}"


def test_no_call_genotype_is_resolved_not_emitted(input_vcf_run):
    """A './.' genotype cannot reach the output: the reads were built from some concrete
    assignment, so the VCF has to report that one."""
    import pysam

    with pysam.VariantFile(str(input_vcf_run.vcf)) as vcf:
        for record in vcf:
            for sample in record.samples.values():
                alleles = sample.get("GT")
                if alleles is not None and all(allele is None for allele in alleles):
                    pytest.fail(f"{record.chrom}:{record.pos} still carries a no-call genotype")


def test_input_variant_reads_are_still_placed_correctly(input_vcf_run):
    """Applying input variants must not disturb where reads say they came from — a large
    insertion in particular used to be unrepresentable in the CIGAR (#326)."""
    failures = invariants.reads_are_placed_where_they_say(input_vcf_run.bam, input_vcf_run.reference)
    assert not failures, (
        f"{len(failures)} misplaced record(s):\n  " + "\n  ".join(failures[:10])
    )


@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_invariants_hold_across_seeds(seed, references, run_simulation):
    """
    The same checks under a different draw order.

    An invariant that only holds for one seed is holding by cancellation. This is the shape of
    #337: the sequence and quality lengths agree in the order the sampler happens to produce, and
    stop agreeing as soon as that order changes. A short-insert configuration is used because it
    has no headroom to absorb an off-by-one.
    """
    run = run_simulation(
        references["uniform"],
        {
            "read_len": 150, "coverage": 8, "paired_ended": True,
            "fragment_mean": 130, "fragment_st_dev": 40,
            "adapters": True, "adapter_preset": "truseq",
            "mutation_rate": 0.005, "rng_seed": seed,
        },
        prefix=f"seed{seed}",
    )
    failures = (
        invariants.reads_are_placed_where_they_say(run.bam, run.reference)
        + invariants.cigars_account_for_every_base(run.bam)
        + invariants.bam_is_sorted_and_indexable(run.bam)
        + invariants.fastq_records_match_the_bam(run.bam, run.fastqs)
    )
    assert not failures, f"[seed {seed}] {len(failures)} failure(s):\n  " + "\n  ".join(failures[:10])
