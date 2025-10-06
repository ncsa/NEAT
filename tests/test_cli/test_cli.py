import argparse
import builtins
import io
import os
from pathlib import Path

import pytest

from neat.cli.cli import Cli, main


def test_cli_registers_read_simulator_subcommand():
    cli = Cli()
    # Argparse stores subparsers in a private map; ensure our command is registered
    subparsers_action = None
    for action in cli.parser._actions:  # type: ignore[attr-defined]
        if isinstance(action, argparse._SubParsersAction):
            subparsers_action = action
            break
    assert subparsers_action is not None, "Subparsers not initialized"
    names = set(subparsers_action._name_parser_map.keys())  # type: ignore[attr-defined]
    assert "read-simulator" in names


def test_main_help_returns_2():
    cli = Cli()
    rc = main(cli.parser, ["--help"])  # argparse SystemExit is caught and returns 2
    assert rc == 2


def test_main_no_subcommand_prints_help_and_returns_1(capsys, tmp_path: Path):
    cli = Cli()
    # Direct log file into a temp directory to avoid leaving logs in repo root
    rc = main(cli.parser, ["--log-name", str(tmp_path / "nolines.log")])
    captured = capsys.readouterr()
    assert rc == 1
    assert "usage:" in captured.out.lower()


def test_logging_creates_named_log_file_and_announces(monkeypatch, tmp_path: Path, capsys):
    cli = Cli()
    logname = tmp_path / "myrun.log"
    # Ensure no file exists beforehand
    if logname.exists():
        logname.unlink()

    # Avoid executing the real runner
    monkeypatch.setattr(
        "neat.cli.commands.read_simulator.read_simulator_runner",
        lambda *args, **kwargs: None,
    )

    rc = main(cli.parser, [
        "--log-name", str(logname),
        # Supply a benign subcommand with minimal required args
        "read-simulator", "-o", str(tmp_path), "-p", "pref"
    ])
    out = capsys.readouterr().out
    # main should create/log the file path and return 0 (success)
    assert f"NEAT run log: {logname.resolve()}" in out
    assert logname.exists()
    assert rc == 0


def test_read_simulator_success_invokes_runner(monkeypatch, tmp_path: Path):
    # Arrange
    cli = Cli()
    called = {}

    def fake_runner(cfg, outdir, prefix):
        called['args'] = (cfg, outdir, prefix)

    # Patch runner used by command
    monkeypatch.setattr("neat.cli.commands.read_simulator.read_simulator_runner", fake_runner)

    cfg = tmp_path / "conf.yml"
    cfg.write_text("reference: ''\n", encoding="utf-8")  # minimal content; not validated here

    rc = main(cli.parser, [
        "--no-log",
        "read-simulator",
        "-c", str(cfg),
        "-o", str(tmp_path),
        "-p", "myprefix",
    ])

    assert rc == 0
    assert called['args'] == (str(cfg), str(tmp_path), "myprefix")


def test_read_simulator_failure_returns_1_and_prints_error(monkeypatch, tmp_path: Path, capsys):
    cli = Cli()

    def boom(*args, **kwargs):
        raise RuntimeError("kaboom")

    monkeypatch.setattr("neat.read_simulator.read_simulator_runner", boom)

    rc = main(cli.parser, [
        "--no-log",
        "read-simulator",
        "-o", str(tmp_path),
        "-p", "x",
    ])

    out = capsys.readouterr().out
    assert rc == 1
    # Error path prints a line starting with 'ERROR:'
    assert "ERROR:" in out
