"""
Description.
"""

import argparse
import gzip
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List
import yaml

__all__ = ["main"]


def print_stderr(msg: str, *, exit_: bool = False):
    print(msg, file=sys.stderr)
    if exit_:
        sys.exit(1)


def gather(paths: Iterable[Path], suffixes: tuple[str, ...]) -> List[Path]:
    out: list[Path] = []
    for p in paths:
        if p.is_file() and p.name.endswith(suffixes):
            out.append(p)
        elif p.is_dir():
            for suf in suffixes:
                out.extend(p.rglob(f"*{suf}"))
    return sorted(set(out))


def is_gzipped(file: Path) -> bool:
    return file.suffix in {".gz", ".bgz"}


def concat(files: List[Path], dest: Path, keep_compression: bool):
    if not files:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)

    if keep_compression:
        if not dest.suffix.endswith("gz"):
            dest = dest.with_suffix(dest.suffix + ".gz")
        with dest.open("wb") as out_f:
            for f in files:
                with f.open("rb") as in_f:
                    shutil.copyfileobj(in_f, out_f)
    else:
        with dest.open("wb") as out_f:
            for f in files:
                opener = gzip.open if is_gzipped(f) else open
                with opener(f, "rb") as in_f:
                    shutil.copyfileobj(in_f, out_f)


def merge_bam(bams: List[Path], dest: Path, samtools: str):
    if not bams:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    samtools = samtools if Path(samtools).exists() else shutil.which(samtools) or samtools
    unsorted = dest.with_suffix(".unsorted.bam")
    cmd_merge = [samtools, "merge", "--no-PG", "-f", str(unsorted), *map(str, bams)]
    subprocess.check_call(cmd_merge)
    subprocess.check_call([samtools, "sort", "-o", str(dest), str(unsorted)])
    subprocess.check_call([samtools, "index", str(dest)])
    unsorted.unlink(missing_ok=True)
    normalize_bam_header(dest, samtools)


def merge_vcf(vcfs: List[Path], dest: Path):
    if not vcfs:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    first, *rest = vcfs
    shutil.copy(first, dest)
    with dest.open("ab") as out_f:
        for vcf in rest:
            with vcf.open("rb") as fh:
                for line in fh:
                    if not line.startswith(b"#"):
                        out_f.write(line)


def normalize_bam_header(bam: Path, samtools: str):
    """Make BAM header deterministic: drop @PG, sort @SQ/@RG."""
    hdr_text = subprocess.check_output([samtools, "view", "-H", str(bam)], text=True)

    hd = []
    sq = []
    rg = []
    other = []
    for line in hdr_text.splitlines():
        if line.startswith("@PG"):
            continue
        if line.startswith("@HD"):
            hd.append(line)
        elif line.startswith("@SQ"):
            sq.append(line)
        elif line.startswith("@RG"):
            # optional: drop volatile CL: field
            parts = [t for t in line.split("\t") if not t.startswith("CL:")]
            line = "\t".join(parts)
            rg.append(line)
        else:
            other.append(line)

    def sn(x):
        for t in x.split("\t"):
            if t.startswith("SN:"):
                return t[3:]
        return x

    def rgid(x):
        for t in x.split("\t"):
            if t.startswith("ID:"):
                return t[3:]
        return x

    sq.sort(key=sn)
    rg.sort(key=rgid)

    new_hdr = "\n".join([*(hd or ["@HD\tVN:1.6\tSO:coordinate"]), *sq, *rg, *other]) + "\n"
    hdr_tmp = bam.with_suffix(".reheader.sam")
    out_tmp = bam.with_suffix(".reheader.bam")
    hdr_tmp.write_text(new_hdr)

    with out_tmp.open("wb") as out_f:
        subprocess.check_call([samtools, "reheader", str(hdr_tmp), str(bam)], stdout=out_f)

    out_tmp.replace(bam)
    hdr_tmp.unlink(missing_ok=True)
    subprocess.check_call([samtools, "index", str(bam)])


def natural_key(path: Path):
    s = path.name
    parts = re.split(r"(\d+)", s)
    return [int(p) if p.isdigit() else p for p in parts]


def order_by_manifest(files: List[Path], manifest_path: Path, suffixes: tuple[str, ...]) -> List[Path]:
    """
    Order files according to entry stems, matched by stem either in the file basename or in its parent directory name.
    """
    m = yaml.safe_load(manifest_path.read_text())
    stems = [Path(e["stem"]).name for e in m.get("entries", [])]

    # Group incoming files by a detected stem token
    by_stem: dict[str, list[Path]] = {}
    for f in files:
        name = f.name
        parent = f.parent.name

        # Try to find a matching stem token in either
        key = None
        for s in stems:
            if name.startswith(s) or parent == s:
                key = s
                break
        if key is None:

            # Strip known suffix and use basename as key
            for suf in suffixes:
                if name.endswith(suf):
                    key = name[: -len(suf)]
                    break
        by_stem.setdefault(key or name, []).append(f)

    ordered: list[Path] = []
    for s in stems:
        ordered.extend(sorted(by_stem.get(s, []), key=natural_key))

    # Include remaining items in a stable order
    remaining = [f for k, v in by_stem.items() if k not in stems for f in v]
    ordered.extend(sorted(remaining, key=natural_key))
    return ordered


def parse_args(argv: list[str] | None = None):
    p = argparse.ArgumentParser(description="Stitch NEAT split-run outputs into one dataset.")
    p.add_argument("-i", "--inputs", nargs="+", type=Path, required=True,
                   help="Input split-run directories or file prefixes")
    p.add_argument("-o", "--output-prefix", type=Path, required=True,
                   help="Prefix for stitched outputs (no extension)")
    p.add_argument("-c", "--config", type=Path, required=True,
                   help="Original NEAT YAML config to detect paired-ended, etc.")
    p.add_argument("--samtools", default="samtools", help="Path to samtools binary")
    p.add_argument("--manifest", type=Path, help="Optional manifest.yaml from split step to enforce ordering")
    return p.parse_args(argv)


def main(argv: list[str] | None = None):
    args = parse_args(argv)

    cfg = yaml.safe_load(args.config.read_text())
    paired = bool(cfg.get("paired_ended", False))

    # These cover both NEAT internal naming (_r1_paired.fq.bgz) and CLI naming (_r1.fastq.gz)
    r1_suffixes = ("_r1_paired.fq", "_r1_paired.fq.gz", "_r1_paired.fq.bgz", "_r1.fastq.gz")
    r2_suffixes = ("_r2_paired.fq", "_r2_paired.fq.gz", "_r2_paired.fq.bgz", "_r2.fastq.gz")
    s1_suffixes = ("_r1_single.fq", "_r1_single.fq.gz", "_r1_single.fq.bgz")  # single-end lane

    fq_r1 = gather(args.inputs, r1_suffixes)
    fq_r2 = gather(args.inputs, r2_suffixes) if paired else []
    fq_single = gather(args.inputs, s1_suffixes)
    bams = gather(args.inputs, ("_golden.bam",))
    vcfs = gather(args.inputs, ("_golden.vcf",))

    if not any((fq_r1, fq_r2, fq_single, bams, vcfs)):
        print_stderr("No NEAT output files found under supplied inputs", exit_=True)

    for lane in (fq_r1, fq_r2):
        if lane:
            gz = [is_gzipped(f) for f in lane]
            if any(gz) and not all(gz):
                print_stderr("Mixing compressed & uncompressed FASTQs; aborting stitching", exit_=True)
    keep_comp = fq_r1 and all(is_gzipped(f) for f in fq_r1)

    # Enforce stable ordering
    if args.manifest and args.manifest.exists():
        fq_r1 = order_by_manifest(fq_r1, args.manifest, r1_suffixes)
        fq_r2 = order_by_manifest(fq_r2, args.manifest, r2_suffixes) if paired else []
        fq_single = order_by_manifest(fq_single, args.manifest, s1_suffixes)
        bams = order_by_manifest(bams, args.manifest, ("_golden.bam",))
        vcfs = order_by_manifest(vcfs, args.manifest, ("_golden.vcf",))
    else:
        fq_r1.sort(key=natural_key)
        fq_r2.sort(key=natural_key)
        fq_single.sort(key=natural_key)
        bams.sort(key=natural_key)
        vcfs.sort(key=natural_key)

    prefix = args.output_prefix
    concat(fq_r1, prefix.with_name(prefix.name + "_read1.fq"), keep_comp)
    if paired:
        concat(fq_r2, prefix.with_name(prefix.name + "_read2.fq"), keep_comp)
    concat(fq_single, prefix.with_name(prefix.name + "_readSingle.fq"), keep_comp)

    merge_bam(bams, prefix.with_name(prefix.name + "_golden.bam"), args.samtools)
    merge_vcf(vcfs, prefix.with_name(prefix.name + "_golden.vcf"))

    print_stderr("Stitching complete!")


if __name__ == "__main__":
    main()
