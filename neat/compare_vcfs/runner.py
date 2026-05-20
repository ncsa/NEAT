"""
Runner for the `neat compare-vcfs` subcommand (issue #297).

Compares a downstream variant caller's VCF against a NEAT-simulated truth VCF
and attributes false negatives to the simulator's own configuration (mutation
bed, target bed, simulated contigs) by reading `simulation_summary.json`.

This module currently implements **only the input-validation scaffold**:
- Validates input file paths.
- Locates `hap.py` on PATH (or via `--happy-bin`).
- Loads and version-checks `simulation_summary.json`.
- Creates the output directory.

The actual hap.py invocation, output parsing, FN attribution, and report
generation are tracked as follow-up work on the same issue.
"""
import json
import logging
import shutil
from pathlib import Path

from ..common import validate_input_path
from ..common.chrom_names import load_chrom_aliases
from ..read_simulator.utils.simulation_summary import SCHEMA_VERSION
from .attribution import attribute_fns, detect_chrom_naming_mismatches
from .happy import run_happy, parse_happy_output
from .reports import (
    build_comparison_summary,
    summarize_fn_reasons,
    write_comparison_summary_json,
    write_comparison_summary_txt,
    write_fn_attribution_plot,
    write_fn_with_reasons,
)

__all__ = ["compare_vcfs_runner", "load_simulation_summary", "discover_happy"]

_LOG = logging.getLogger(__name__)

# An install hint surfaced when hap.py is not on PATH and --happy-bin is not given.
_HAPPY_INSTALL_HINT = (
    "hap.py was not found on $PATH. Install via "
    "`conda install -c bioconda hap.py`, "
    "or use the official Docker image, or pass --happy-bin /path/to/hap.py."
)


class HappyNotFoundError(RuntimeError):
    """Raised when hap.py cannot be located."""


class SimulationSummaryError(RuntimeError):
    """Raised when simulation_summary.json is missing, malformed, or incompatible."""


def compare_vcfs_runner(
    golden_vcf: str,
    called_vcf: str,
    neat_run_dir: str,
    output_dir: str,
    reference: str | None = None,
    target_bed: str | None = None,
    happy_bin: str | None = None,
    plot: bool = False,
    chrom_aliases: str | None = None,
):
    """
    Run the comparison pipeline.

    :param golden_vcf: NEAT-simulated truth VCF.
    :param called_vcf: Downstream variant caller's VCF.
    :param neat_run_dir: Directory containing the NEAT simulator output, including
        `simulation_summary.json`.
    :param output_dir: Where to write `comparison_summary.{json,txt}` and
        `FN_with_reasons.vcf`. Created if it does not exist.
    :param reference: Optional path to the reference FASTA (forwarded to hap.py).
    :param target_bed: Optional BED of regions to restrict comparison to.
    :param happy_bin: Optional explicit path to the hap.py binary.
    """
    golden_path = Path(golden_vcf).resolve()
    called_path = Path(called_vcf).resolve()
    run_dir_path = Path(neat_run_dir).resolve()
    out_dir_path = Path(output_dir).resolve()

    validate_input_path(golden_path)
    validate_input_path(called_path)
    if not run_dir_path.is_dir():
        raise FileNotFoundError(f"--neat-run-dir does not exist or is not a directory: {run_dir_path}")

    if reference is not None:
        validate_input_path(Path(reference).resolve())
    if target_bed is not None:
        validate_input_path(Path(target_bed).resolve())
    if chrom_aliases is not None:
        validate_input_path(Path(chrom_aliases).resolve())
    aliases = load_chrom_aliases(chrom_aliases)
    if aliases:
        _LOG.info(f"Loaded {len(aliases)} chrom alias mapping(s) from {chrom_aliases}")

    happy = discover_happy(happy_bin)
    _LOG.info(f"Using hap.py at: {happy}")

    summary = load_simulation_summary(run_dir_path)
    _LOG.info(
        f"Loaded simulation_summary.json (schema v{summary['schema_version']}, "
        f"neat {summary['neat_version']}, "
        f"{summary['delivered'].get('total_variants')} simulated variants "
        f"across {len(summary['delivered'].get('contigs_simulated', []))} contigs)."
    )

    out_dir_path.mkdir(parents=True, exist_ok=True)
    _LOG.info(f"Output directory: {out_dir_path}")

    happy_prefix = out_dir_path / "happy"
    happy_vcf = run_happy(
        happy_bin=happy,
        golden_vcf=golden_path,
        called_vcf=called_path,
        output_prefix=happy_prefix,
        reference=Path(reference).resolve() if reference else None,
        target_bed=Path(target_bed).resolve() if target_bed else None,
    )
    buckets = parse_happy_output(happy_vcf)
    _LOG.info(
        f"hap.py classification: TP={len(buckets['TP'])} "
        f"FN={len(buckets['FN'])} FP={len(buckets['FP'])}"
    )

    chrom_warnings = detect_chrom_naming_mismatches(summary, aliases=aliases)
    for w in chrom_warnings:
        _LOG.warning(w["message"])

    # A fully mismatched BED is unusable for attribution; skip it so we don't
    # report misleading 'outside_<bed>' counts for chroms NEAT never actually
    # checked against the BED. The warning surfaces the underlying cause.
    unusable_beds = {w["bed"] for w in chrom_warnings if w.get("type") == "chrom_naming_mismatch"}

    fn_reasons = attribute_fns(buckets["FN"], summary, aliases=aliases, skip_beds=unusable_beds)
    reason_counts = summarize_fn_reasons(fn_reasons)
    if fn_reasons:
        _LOG.info(f"FN attribution: {reason_counts}")

    fn_with_reasons_path = out_dir_path / "FN_with_reasons.vcf"
    comparison_json_path = out_dir_path / "comparison_summary.json"
    comparison_txt_path = out_dir_path / "comparison_summary.txt"

    write_fn_with_reasons(happy_vcf, fn_reasons, fn_with_reasons_path)
    _LOG.info(f"Wrote {fn_with_reasons_path}")

    counts = {k: len(v) for k, v in buckets.items()}
    report = build_comparison_summary(
        golden_vcf=golden_path,
        called_vcf=called_path,
        neat_run_dir=run_dir_path,
        simulation_summary_path=run_dir_path / "simulation_summary.json",
        happy_output_vcf=happy_vcf,
        happy_output_prefix=happy_prefix,
        counts=counts,
        fn_attribution=reason_counts,
        fn_with_reasons_vcf=fn_with_reasons_path,
        comparison_summary_json=comparison_json_path,
        comparison_summary_txt=comparison_txt_path,
        warnings=chrom_warnings,
    )
    write_comparison_summary_json(report, comparison_json_path)
    write_comparison_summary_txt(report, comparison_txt_path)
    _LOG.info(f"Wrote {comparison_json_path}")
    _LOG.info(f"Wrote {comparison_txt_path}")

    if plot:
        plot_path = out_dir_path / "fn_attribution.png"
        write_fn_attribution_plot(reason_counts, plot_path)
        _LOG.info(f"Wrote {plot_path}")


def discover_happy(explicit_path: str | None) -> Path:
    """
    Resolve the hap.py binary. Explicit path wins; otherwise look on $PATH.

    :param explicit_path: Path passed via --happy-bin, or None.
    :return: Absolute path to a hap.py executable.
    :raises HappyNotFoundError: if neither location yields an executable.
    """
    if explicit_path is not None:
        p = Path(explicit_path).resolve()
        if not p.is_file():
            raise HappyNotFoundError(f"--happy-bin path does not exist: {p}")
        return p
    found = shutil.which("hap.py")
    if found is None:
        raise HappyNotFoundError(_HAPPY_INSTALL_HINT)
    return Path(found).resolve()


def load_simulation_summary(neat_run_dir: Path) -> dict:
    """
    Read and validate `simulation_summary.json` from the NEAT run directory.

    :raises SimulationSummaryError: if the file is missing, malformed, or its
        schema version is incompatible with this build.
    """
    summary_path = neat_run_dir / "simulation_summary.json"
    if not summary_path.is_file():
        raise SimulationSummaryError(
            f"{summary_path} not found. Did the simulator run finish? "
            "Re-run `neat read-simulator` to produce one."
        )
    try:
        with open(summary_path) as fh:
            data = json.load(fh)
    except json.JSONDecodeError as exc:
        raise SimulationSummaryError(f"{summary_path} is not valid JSON: {exc}") from exc

    actual = data.get("schema_version")
    if actual != SCHEMA_VERSION:
        raise SimulationSummaryError(
            f"{summary_path} schema_version is {actual!r}, expected {SCHEMA_VERSION!r}. "
            "Regenerate by re-running `neat read-simulator` with a current NEAT build."
        )
    return data
