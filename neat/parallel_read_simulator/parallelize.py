"""
Split a reference genome, run NEAT simulations in
parallel, and stitch their outputs back together.
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple

from .split_inputs import main as split_main
from .stitch_outputs import main as stitch_main

import logging

__all__ = ["parallelize_main"]

_LOG = logging.getLogger(__name__)


# Parallelization worker
def worker(args: Tuple[List[str], Path]):
    """
    Run one NEAT simulation in its own working directory.

    Each worker is executed in a separate process by the
    :class:`concurrent.futures.ProcessPoolExecutor`.
    """
    cmd, workdir = args
    workdir.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.check_call(cmd, cwd=workdir)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"NEAT job failed (exit {e.returncode}) in {workdir.name}: {' '.join(cmd)}"
        ) from e
    return workdir


# Argument parser
def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    """
    Parse command‑line arguments for the parallel pipeline.
    """
    p = argparse.ArgumentParser(
        description="Split, simulate in parallel, and stitch NEAT outputs back together.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "config",
        type=Path,
        help="NEAT YAML/YML config containing the 'reference:' field",
    )
    p.add_argument(
        "--outdir",
        type=Path,
        required=False,
        default=None,
        help="Top-level directory for splits and stitched results (optional)",
    )

    # Splitting options
    split = p.add_argument_group("splitting options")
    split.add_argument(
        "--by",
        choices=["contig", "size"],
        default=None,
        help="Split mode",
    )
    split.add_argument(
        "--size",
        type=int,
        default=None,
        help="Target chunk size when --by size",
    )
    split.add_argument(
        "--cleanup-splits",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Delete the 'splits' directory after stitching completes",
    )
    split.add_argument(
        "--reuse-splits",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Skip the splitting step and reuse any YAML/FASTA files already present",
    )

    # Simulation options
    sim = p.add_argument_group("simulation options")
    sim.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="Maximum number of parallel NEAT jobs to run",
    )
    sim.add_argument(
        "--neat-cmd",
        default=None,
        help="Command to launch the read simulator",
    )

    # Stitching options
    stitch = p.add_argument_group("stitching options")
    stitch.add_argument(
        "--samtools",
        default=None,
        help="Path to samtools executable used by stitch_outputs.py",
    )
    stitch.add_argument(
        "--final-prefix",
        type=Path,
        default=None,
        help="Prefix for stitched outputs",
    )

    return p.parse_args(argv)


def parallelize_main(argv: List[str] | None = None) -> None:
    """
    This function executes the three high-level steps of parallelization:

    1. Split the reference sequence into chunks
    2. Launch NEAT's read simulator on each chunk
    3. Stitch the resulting outputs back into a single dataset
    """
    args = parse_args(argv)
    try:
        import yaml
    except ModuleNotFoundError:
        _LOG.error("PyYAML not installed — cannot read config")
        sys.exit(1)

    # Load the YAML configuration file
    cfg_dict = yaml.safe_load(args.config.read_text())

    # Helper functions: treat '.' as unset and cast values
    sentinel = {None, "."}

    def yaml_val(key: str):
        v = cfg_dict.get(key)
        return None if v in sentinel else v

    def cfg_val(key: str, typ: type | None = None):
        v = yaml_val(key)
        if v is None:
            return None
        if typ is Path:
            return Path(str(v))
        if typ and not isinstance(v, typ):
            return typ(v)
        return v

    def coalesce(cli_val, key: str, typ: type | None, default):
        if cli_val is not None:
            return cli_val
        y = cfg_val(key, typ)
        return y if y is not None else default

    def coalesce_bool(cli_val, key: str, default: bool = False) -> bool:
        if cli_val is not None:
            return bool(cli_val)
        y = yaml_val(key)
        if isinstance(y, bool):
            return y
        if isinstance(y, str):
            s = y.strip().lower()
            if s in {"1", "true", "t", "yes", "y", "on"}:
                return True
            if s in {"0", "false", "f", "no", "n", "off"}:
                return False
        return default

    # Resolve reference path issues early

    ref_path = cfg_val("reference", Path)
    if not ref_path:
        _LOG.error("'reference:' missing in config")
        sys.exit(1)

    ref_path = ref_path.expanduser().resolve()
    if not ref_path.is_file():
        _LOG.error(f"'reference:' not found or file does not exist: {ref_path}")
        sys.exit(1)

    # Merge config values

    # Splitting
    args.by = coalesce(args.by, "by", str, "contig")
    args.size = coalesce(args.size, "size", int, 500000 if args.by == "size" else None)

    # Simulate
    args.jobs = coalesce(args.jobs, "jobs", int, os.cpu_count() or 2)
    args.neat_cmd = coalesce(args.neat_cmd, "neat_cmd", str, "neat read-simulator")

    # Stitching
    args.samtools = coalesce(args.samtools, "samtools", str, "samtools")
    args.final_prefix = coalesce(args.final_prefix, "final_prefix", Path, Path("stitched/final"))

    # Outputs
    cwd = Path.cwd()
    default_outdir = (cwd / f"{Path(args.config).stem}_parallel").resolve()
    args.outdir = coalesce(args.outdir, "outdir", Path, default_outdir)
    if not args.outdir.is_absolute():
        args.outdir = (cwd / args.outdir).resolve()

    # Boolean toggles
    args.cleanup_splits = coalesce_bool(args.cleanup_splits, "cleanup_splits", False)
    args.reuse_splits = coalesce_bool(args.reuse_splits, "reuse_splits", False)

    # Normalize final_prefix to absolute
    if not args.final_prefix.is_absolute():
        args.final_prefix = (args.outdir / args.final_prefix).resolve()

    # Create directories after resolution
    args.outdir.mkdir(parents=True, exist_ok=True)
    splits_dir = args.outdir / "splits"
    sims_dir = args.outdir / "sim_runs"
    splits_dir.mkdir(parents=True, exist_ok=True)
    sims_dir.mkdir(parents=True, exist_ok=True)

    # 1) Split (or reuse)
    need_split = not args.reuse_splits

    if args.reuse_splits and list(splits_dir.glob("*.y*ml")):
        _LOG.info("[parallel] Recycling existing split FASTA/YAML files – skipping split step")
        need_split = False

    if need_split:

        t0 = time.time()
        _LOG.info("[parallel] Splitting reference...")
        split_argv = [
            str(ref_path),
            str(args.config),
            "--outdir", str(splits_dir),
            "--by", args.by,
        ]

        if args.by == "size":
            split_argv += ["--size", str(args.size)]
        split_main(split_argv)
        split_sec = time.time() - t0

    else:
        split_sec = 0.0

    # 2) Submit NEAT sims in a deterministic order using the manifest
    manifest_path = splits_dir / "manifest.yaml"
    if manifest_path.exists():
        m = yaml.safe_load(manifest_path.read_text())
        yaml_files = [Path(e["yaml"]) for e in m.get("entries", [])]
    else:
        yaml_files = sorted(list(splits_dir.glob("*.yaml")) + list(splits_dir.glob("*.yml")))

    if not yaml_files:
        _LOG.error("No split YAML/YML configs found — aborting")
        sys.exit(1)

    _LOG.info(f"[parallel] Launching {len(yaml_files)} NEAT job(s) (max {args.jobs} in parallel)...")
    t1 = time.time()
    futures = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for yaml_path in yaml_files:
            stem = yaml_path.stem
            workdir = sims_dir / stem
            neat_cmd = args.neat_cmd.split()
            # Launch NEAT via CLI (default 'neat read-simulator')
            cmd = [*neat_cmd, "-c", str(yaml_path.resolve()), "-o", stem]
            futures.append(pool.submit(worker, (cmd, workdir)))

        for fut in as_completed(futures):
            try:
                wd = fut.result()
                _LOG.info(f" Finished {wd.name}")
            except Exception as exc:
                _LOG.error(f" Simulation failed: {exc}")
                sys.exit(1)
    sim_sec = time.time() - t1

    # 3) Stitch (also in deterministic order, via --manifest when available)

    _LOG.info("[parallel] Stitching per-chunk outputs...")
    t2 = time.time()

    stitch_argv = [
        "-i", str(sims_dir),
        "-o", str(args.final_prefix),
        "-c", str(args.config),
        "--samtools", args.samtools,
    ]

    if manifest_path.exists():
        stitch_argv += ["--manifest", str(manifest_path)]

    stitch_main(stitch_argv)
    stitch_sec = time.time() - t2

    if args.cleanup_splits:
        _LOG.info("[parallel] Cleaning up split FASTA/configs...")
        shutil.rmtree(splits_dir, ignore_errors=True)

    _LOG.info(
        f"[parallel] Pipeline complete — stitched files under {args.final_prefix.parent}\n"
        f"Timings:\n"
        f"  split : {split_sec:6.1f} s\n"
        f"  sim   : {sim_sec:6.1f} s\n"
        f"  stitch: {stitch_sec:6.1f} s\n"
        f"  total : {split_sec + sim_sec + stitch_sec:6.1f} s"
    )
