"""
CLI-layer tests for `neat compare-vcfs` — subcommand registration, argument
parsing, and routing into the runner.
"""
import argparse

import pytest

from neat.cli.cli import Cli
from neat.cli.commands.compare_vcfs import Command


def _compare_vcfs_subparser():
    cli = Cli()
    subparsers_action = next(
        a for a in cli.parser._actions if isinstance(a, argparse._SubParsersAction)
    )
    return subparsers_action._name_parser_map.get("compare-vcfs")


# ===========================================================================
# Subcommand registration
# ===========================================================================

def test_compare_vcfs_subcommand_is_registered():
    """`neat compare-vcfs` must be discovered by Cli's pkgutil-based scan."""
    sub = _compare_vcfs_subparser()
    assert sub is not None, "compare-vcfs subcommand is not registered"


def test_compare_vcfs_help_text_lists_all_documented_flags():
    """Regression guard: every flag described in the docs must appear in --help."""
    sub = _compare_vcfs_subparser()
    help_text = sub.format_help()
    for flag in [
        "--neat-run-dir",
        "--output-dir",
        "--reference",
        "--target-bed",
        "--happy-bin",
        "--plot",
        "--chrom-aliases",
    ]:
        assert flag in help_text, f"--help missing {flag}"


# ===========================================================================
# Argument routing — Command.execute → compare_vcfs_runner
# ===========================================================================

def test_command_routes_all_args_to_runner(monkeypatch):
    """Parsed args must reach compare_vcfs_runner with the right kwargs."""
    captured = {}

    def fake_runner(**kwargs):
        captured.update(kwargs)

    from neat.cli.commands import compare_vcfs as cmd_mod
    monkeypatch.setattr(cmd_mod, "compare_vcfs_runner", fake_runner)

    parser = argparse.ArgumentParser()
    cmd = Command(parser)
    args = parser.parse_args([
        "/g.vcf", "/c.vcf",
        "--neat-run-dir", "/run",
        "--output-dir", "/out",
        "--reference", "/ref.fa",
        "--target-bed", "/t.bed",
        "--happy-bin", "/bin/hap.py",
        "--plot",
        "--chrom-aliases", "/aliases.tsv",
    ])
    cmd.execute(args)

    assert captured == {
        "golden_vcf": "/g.vcf",
        "called_vcf": "/c.vcf",
        "neat_run_dir": "/run",
        "output_dir": "/out",
        "reference": "/ref.fa",
        "target_bed": "/t.bed",
        "happy_bin": "/bin/hap.py",
        "plot": True,
        "chrom_aliases": "/aliases.tsv",
    }


def test_command_optional_flags_default_to_none_or_false():
    """When optional flags aren't supplied, parsed values match runner defaults."""
    parser = argparse.ArgumentParser()
    Command(parser)
    args = parser.parse_args([
        "/g.vcf", "/c.vcf",
        "--neat-run-dir", "/run",
        "--output-dir", "/out",
    ])
    assert args.reference is None
    assert args.target_bed is None
    assert args.happy_bin is None
    assert args.plot is False
    assert args.chrom_aliases is None
