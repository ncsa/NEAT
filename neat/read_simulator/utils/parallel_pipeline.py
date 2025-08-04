from __future__ import annotations

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
    p.add_argument("--outdir", type=Path, required=True, help="Top‑level directory to hold splits, per‑chunk outputs, and final stitched results",)

    split = p.add_argument_group("splitting options")
    split.add_argument("--by", choices=["contig", "size"], default="contig", help="Split mode")
    split.add_argument("--size", type=int, default=1_000_000, help="Target chunk size when --by size")
    split.add_argument("--cleanup-splits", action="store_true", help="Delete the 'splits' directory after stitching completes")
    split.add_argument("--reuse-splits", action="store_true", help="Skip the splitting step and reuse any YAML/FASTA files already present in the 'splits' folder")

    sim = p.add_argument_group("simulation options")
    sim.add_argument("--jobs", type=int, default=os.cpu_count() or 2, help="Maximum number of parallel NEAT jobs to run",)
    sim.add_argument("--neat-cmd", default="python -m neat.read_simulator.run", help="Command used to launch the read simulator (will be split on whitespace)",)

    stitch = p.add_argument_group("stitching options")
    stitch.add_argument("--samtools", default="samtools", help="Path to samtools executable used by stitch_outputs.py",)
    stitch.add_argument("--final-prefix", type=Path, default=Path("stitched/final"), help="Prefix (no extension) for the stitched outputs",)

    return p.parse_args(argv)


def main(argv: List[str] | None = None):
    args = parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    splits_dir = args.outdir / "splits"
    sims_dir = args.outdir / "sim_runs"
    splits_dir.mkdir(parents=True, exist_ok=True)
    sims_dir.mkdir(parents=True, exist_ok=True)

    final_prefix = args.final_prefix if args.final_prefix.is_absolute() else args.outdir / args.final_prefix

    try:
        import yaml
    except ModuleNotFoundError:
        print("PyYAML not installed — cannot read config", file=sys.stderr)
        sys.exit(1)

    cfg_dict = yaml.safe_load(args.config.read_text())
    ref_path = Path(cfg_dict.get("reference", "")).expanduser().resolve()

    if not ref_path.is_file():
        print(f"'reference:' not found or file does not exist: {ref_path}", file=sys.stderr)
        sys.exit(1)

    # 1. Split the input reference

    need_split = not args.reuse_splits

    if args.reuse_splits and list((args.outdir / "splits").glob("*.y*ml")):
        print("[parallel] Recycling existing split FASTA/YAML files – skipping split step")
    else:
        need_split = True

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

    # 2. Run NEAT simulations in parallel

    yaml_files = sorted(list(splits_dir.glob("*.yaml")) + list(splits_dir.glob("*.yml")))

    if not yaml_files:
        print("No split YAML/YML configs generated — aborting")
        sys.exit(1)

    print(f"[parallel] Launching {len(yaml_files)} NEAT job(s) (max {args.jobs} in parallel)...")
    futures = []

    t1 = time.time()
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for yaml_path in yaml_files:
            stem = yaml_path.stem
            workdir = sims_dir / stem
            neat_cmd = args.neat_cmd.split()
            cmd = [*neat_cmd, "-c", str(yaml_path.resolve()), "-o", str(workdir)]
            futures.append(pool.submit(worker, (cmd, workdir)))

        # progress / error handling
        for fut in as_completed(futures):
            try:
                wd = fut.result()
                print(f" Finished {wd.name}")
            except Exception as exc:
                print(f" Simulation failed: {exc}")
                sys.exit(1)

    sim_sec = time.time() - t1

    # 3. Stitch outputs
    stitch_cmd = [
        sys.executable,
        str(Path(__file__).with_name("stitch_outputs.py")),
        "-i", *map(str, sims_dir.iterdir()),
        "-o", str(final_prefix),
        "-c", str(args.config),
        "--samtools", args.samtools,
    ]

    print("[parallel] Stitching per‑chunk outputs...")
    t2 = time.time()

    try:
        subprocess.check_call(stitch_cmd)
    except subprocess.CalledProcessError as e:
        print(f"Stitch step failed (exit {e.returncode})")
        sys.exit(1)

    print(f"[parallel] Pipeline complete — stitched files under {final_prefix.parent}")

    stitch_sec = time.time() - t2
    total_sec = time.time() - t0

    if args.cleanup_splits:
        print("[parallel] Cleaning up split FASTA/configs...")
        shutil.rmtree(splits_dir, ignore_errors=True)

    print(f"[parallel] Pipeline complete — stitched files under {final_prefix.parent}\n"
          f"Timings:\n"
          f"  split : {split_sec:6.1f} s\n"
          f"  sim   : {sim_sec:6.1f} s\n"
          f"  stitch: {stitch_sec:6.1f} s\n"
          f"  total : {total_sec:6.1f} s"
          )


if __name__ == "__main__":
    main()
