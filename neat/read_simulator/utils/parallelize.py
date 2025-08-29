"""
Description.
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List


# Parallelization worker
def worker(args: tuple[List[str], Path]):
    """Run one NEAT simulation in its own working directory."""
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
def parse_args(argv: List[str] | None = None):
    p = argparse.ArgumentParser(
        description="Split, simulate in parallel, and stitch NEAT outputs back together.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("config", type=Path, help="NEAT YAML/YML config containing the 'reference:' field")
    p.add_argument("--outdir", type=Path, required=False, default=None,
                   help="Top-level directory for splits and stitched results (optional)")

    split = p.add_argument_group("splitting options")
    split.add_argument("--by", choices=["contig", "size"], default=None, help="Split mode")
    split.add_argument("--size", type=int, default=None, help="Target chunk size when --by size")
    split.add_argument("--cleanup-splits", action=argparse.BooleanOptionalAction, default=None,
                       help="Delete the 'splits' directory after stitching completes")
    split.add_argument("--reuse-splits", action=argparse.BooleanOptionalAction, default=None,
                       help="Skip the splitting step and reuse any YAML/FASTA files already present")

    sim = p.add_argument_group("simulation options")
    sim.add_argument("--jobs", type=int, default=None, help="Maximum number of parallel NEAT jobs to run")
    sim.add_argument("--neat-cmd", default=None, help="Command to launch the read simulator")

    stitch = p.add_argument_group("stitching options")
    stitch.add_argument("--samtools", default=None, help="Path to samtools executable used by stitch_outputs.py")
    stitch.add_argument("--final-prefix", type=Path, default=None, help="Prefix for stitched outputs")

    return p.parse_args(argv)


def main(argv: List[str] | None = None):

    args = parse_args(argv)
    try:
        import yaml
    except ModuleNotFoundError:
        print("PyYAML not installed — cannot read config", file=sys.stderr)
        sys.exit(1)

    cfg_dict = yaml.safe_load(args.config.read_text())

    # Helper functions: treat '.' as unset and cast values
    sentinel = {None, "."}

    def yaml_val(key):
        v = cfg_dict.get(key)
        return None if v in sentinel else v

    def cfg_val(key, typ=None):
        v = yaml_val(key)
        if v is None:
            return None
        if typ is Path:
            return Path(str(v))
        if typ and not isinstance(v, typ):
            return typ(v)
        return v

    def coalesce(cli_val, key, typ, default):
        if cli_val is not None:
            return cli_val
        y = cfg_val(key, typ)
        return y if y is not None else default

    def coalesce_bool(cli_val, key, default=False):
        if cli_val is not None:
            return cli_val
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
        print("'reference:' missing in config", file=sys.stderr)
        sys.exit(1)

    ref_path = ref_path.expanduser().resolve()
    if not ref_path.is_file():
        print(f"'reference:' not found or file does not exist: {ref_path}", file=sys.stderr)
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
        print("[parallel] Recycling existing split FASTA/YAML files – skipping split step")
        need_split = False

    if need_split:
        split_cmd = [
            sys.executable,
            str(Path(__file__).with_name("split_inputs.py")),
            str(ref_path),
            str(args.config),
            "--outdir", str(splits_dir),
            "--by", args.by,
        ]
        if args.by == "size":
            split_cmd += ["--size", str(args.size)]

        t0 = time.time()
        print("[parallel] Splitting reference...")
        subprocess.check_call(split_cmd)
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
        print("No split YAML/YML configs found — aborting")
        sys.exit(1)

    print(f"[parallel] Launching {len(yaml_files)} NEAT job(s) (max {args.jobs} in parallel)...")
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
                print(f" Finished {wd.name}")
            except Exception as exc:
                print(f" Simulation failed: {exc}")
                sys.exit(1)

    sim_sec = time.time() - t1

    # 3) Stitch (also in deterministic order, via --manifest when available)
    stitch_cmd = [
        sys.executable,
        str(Path(__file__).with_name("stitch_outputs.py")),
        "-i", str(sims_dir),
        "-o", str(args.final_prefix),
        "-c", str(args.config),
        "--samtools", args.samtools,
    ]
    if manifest_path.exists():
        stitch_cmd += ["--manifest", str(manifest_path)]

    print("[parallel] Stitching per-chunk outputs...")
    t2 = time.time()
    try:
        subprocess.check_call(stitch_cmd)
    except subprocess.CalledProcessError as e:
        print(f"Stitch step failed (exit {e.returncode})")
        sys.exit(1)

    stitch_sec = time.time() - t2

    if args.cleanup_splits:
        print("[parallel] Cleaning up split FASTA/configs...")
        shutil.rmtree(splits_dir, ignore_errors=True)

    print(f"[parallel] Pipeline complete — stitched files under {args.final_prefix.parent}\n"
          f"Timings:\n"
          f"  split : {split_sec:6.1f} s\n"
          f"  sim   : {sim_sec:6.1f} s\n"
          f"  stitch: {stitch_sec:6.1f} s\n"
          f"  total : {split_sec + sim_sec + stitch_sec:6.1f} s")


if __name__ == "__main__":
    main()
