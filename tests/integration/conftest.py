"""
Fixtures for the integration suite: generated references, and a helper that runs the simulator
the way a user does.

These tests are skipped unless NEAT_INTEGRATION is set, because each one runs a real simulation
and the matrix takes minutes rather than the second and a half the unit suite takes. See
tests/integration/README.md.
"""

from __future__ import annotations

import os
import random
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]

_SKIP_REASON = (
    "integration suite: set NEAT_INTEGRATION=1 to run (each test runs a full simulation)"
)

requires_integration = pytest.mark.skipif(
    not os.environ.get("NEAT_INTEGRATION"), reason=_SKIP_REASON
)


def _write_fasta(path: Path, contigs: list[tuple[str, str]]) -> Path:
    with path.open("w") as handle:
        for name, sequence in contigs:
            handle.write(f">{name}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i:i + 60] + "\n")
    return path


@pytest.fixture(scope="session")
def references(tmp_path_factory) -> dict[str, Path]:
    """
    Generated rather than committed, so the suite carries no large fixture files and the
    composition of each reference is visible here.

    Length is a tradeoff: long enough that a defect affecting a fraction of a percent of reads
    shows up at all, short enough that a config finishes in seconds.
    """
    directory = tmp_path_factory.mktemp("references")
    rng = random.Random(20260828)

    uniform = "".join(rng.choice("ACGT") for _ in range(120_000))

    gc_rng = random.Random(11)
    gc_rich = "".join(gc_rng.choices("GCGCGCAT", weights=[3, 3, 3, 3, 1, 1, 1, 1], k=120_000))

    rep_rng = random.Random(13)
    blocks = []
    for _ in range(300):
        blocks.append("".join(rep_rng.choice("ACGT") for _ in range(200)))
        blocks.append("CAGGTA" * 50)
    repetitive = "".join(blocks)[:120_000]

    multi_rng = random.Random(14)
    multi = [
        (f"ctg{i}", "".join(multi_rng.choice("ACGT") for _ in range(length)))
        for i, length in enumerate((60_000, 30_000, 18_000, 12_000))
    ]

    return {
        "uniform": _write_fasta(directory / "uniform.fa", [("chr1", uniform)]),
        "gc_rich": _write_fasta(directory / "gc_rich.fa", [("chr1", gc_rich)]),
        "repetitive": _write_fasta(directory / "repetitive.fa", [("chr1", repetitive)]),
        "multi_contig": _write_fasta(directory / "multi.fa", multi),
    }


@pytest.fixture(scope="session")
def mixed_variant_vcf(tmp_path_factory, references) -> Path:
    """An input VCF covering every branch of parse_input_vcf: SNV, insertion, deletion, an MNP
    (which becomes an UnknownVariant), a no-call genotype, and a multiallelic record. Declares
    its own INFO and FILTER keys, and uses one it does not declare, since both paths have to
    survive into the golden VCF."""
    from Bio import SeqIO

    reference = str(next(SeqIO.parse(str(references["uniform"]), "fasta")).seq).upper()
    complement = str.maketrans("ACGT", "TGCA")
    path = tmp_path_factory.mktemp("input_vcf") / "mixed.vcf"

    def at(position):                       # 1-based VCF position -> 0-based reference index
        return reference[position - 1]

    rows = [
        f"chr1\t2000\trs_snv\t{at(2000)}\t{'A' if at(2000) != 'A' else 'C'}\t50\tPASS\tAF=0.5\tGT:DP\t0|1:30",
        f"chr1\t8000\trs_nocall\t{at(8000)}\t{'A' if at(8000) != 'A' else 'C'}\t50\tPASS\tAF=0.5\tGT:DP\t./.:22",
        f"chr1\t14000\trs_del\t{reference[13999:14005]}\t{at(14000)}\t60\tLowQual\tAF=0.25\tGT:DP\t1|1:18",
        f"chr1\t20000\trs_ins\t{at(20000)}\t{at(20000)}GGGGTTTT\t70\tPASS\tDB\tGT:DP\t1|0:25",
        f"chr1\t26000\trs_mnp\t{reference[25999:26002]}\t"
        f"{reference[25999:26002].translate(complement)}\t80\tPASS\tAF=1.0\tGT:DP\t1|1:27",
        f"chr1\t32000\trs_multi\t{at(32000)}\t"
        f"{'A' if at(32000) != 'A' else 'C'},{'T' if at(32000) != 'T' else 'G'}\t90\tPASS\tAF=0.3\tGT:DP\t1|2:31",
    ]
    path.write_text(
        "##fileformat=VCFv4.2\n"
        '##FILTER=<ID=LowQual,Description="Low quality">\n'
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
        # DB is deliberately left undeclared: NEAT has to synthesise a declaration for it.
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        + "\n".join(rows) + "\n"
    )
    return path


class SimulationRun:
    """The output of one simulation, with the paths the invariant checks need."""

    def __init__(self, output_dir: Path, reference: Path, prefix: str = "sim"):
        self.output_dir = output_dir
        self.reference = reference
        self.prefix = prefix

    @property
    def bam(self) -> Path:
        return self.output_dir / f"{self.prefix}_golden.bam"

    @property
    def vcf(self) -> Path:
        return self.output_dir / f"{self.prefix}_golden.vcf.gz"

    @property
    def fastqs(self) -> list[Path]:
        return sorted(self.output_dir.glob(f"{self.prefix}*.fastq.gz"))


@pytest.fixture(scope="session")
def run_simulation(tmp_path_factory):
    """
    Run the simulator in a subprocess, the way a user invokes it, and fail loudly if it does not
    exit 0 — a run that dies partway is itself a finding, and several of the defects this suite
    exists for did exactly that after writing the FASTQs.

    PYTHONPATH is set to the repo so this works against a working tree that has not been pip
    installed, which is also why tests/test_cli/test_basic_cli.py fails locally for most people.
    """
    counter = {"n": 0}

    def _run(reference: Path, config: dict, prefix: str = "sim") -> SimulationRun:
        counter["n"] += 1
        directory = tmp_path_factory.mktemp(f"run{counter['n']}")
        output_dir = directory / "out"
        output_dir.mkdir()

        settings = {
            "reference": str(reference),
            "ploidy": 2,
            "produce_bam": True,
            "produce_fastq": True,
            "produce_vcf": False,
            "rng_seed": 20260828,
            "threads": 1,
            "overwrite_output": True,
        }
        settings.update(config)

        config_path = directory / "config.yml"
        config_path.write_text(yaml.safe_dump(settings, sort_keys=False))

        env = dict(os.environ, PYTHONPATH=str(REPO_ROOT), PYTHONHASHSEED="0")
        result = subprocess.run(
            [sys.executable, "-m", "neat", "--no-log", "read-simulator",
             "-c", str(config_path), "-o", str(output_dir), "-p", prefix],
            capture_output=True, text=True, env=env, cwd=str(directory),
        )
        if result.returncode != 0:
            tail = "\n".join((result.stderr or result.stdout).splitlines()[-25:])
            raise AssertionError(
                f"simulation exited {result.returncode} for config {config}\n{tail}"
            )
        return SimulationRun(output_dir, reference, prefix)

    return _run
