"""
Tests for neat/quality_score_modeling/presets.py and the --quality-preset
CLI wiring in neat/cli/commands/model_qual_score.py.
"""
import argparse

import pytest

from neat.quality_score_modeling.presets import QUALITY_PRESETS


# ===========================================================================
# Preset data integrity
# ===========================================================================

def test_all_presets_are_sorted():
    for name, bins in QUALITY_PRESETS.items():
        assert bins == sorted(bins), f"preset {name!r} bins are not sorted"


def test_all_preset_bins_are_positive():
    for name, bins in QUALITY_PRESETS.items():
        assert all(b > 0 for b in bins), f"preset {name!r} has non-positive bin"


def test_novaseq_has_four_bins():
    assert len(QUALITY_PRESETS["novaseq"]) == 4


def test_nextseq2000_has_four_bins():
    assert len(QUALITY_PRESETS["nextseq2000"]) == 4


def test_nextseq500_has_five_bins():
    assert len(QUALITY_PRESETS["nextseq500"]) == 5


def test_novaseq_bins_match_illumina_spec():
    assert QUALITY_PRESETS["novaseq"] == [2, 12, 23, 37]


def test_nextseq2000_bins_match_illumina_spec():
    assert QUALITY_PRESETS["nextseq2000"] == [2, 12, 26, 37]


# ===========================================================================
# CLI --quality-preset wiring
# ===========================================================================

def _parse_qual_score_args(argv):
    """Run the model-qual-score argument parser against argv."""
    from neat.cli.commands.model_qual_score import Command
    parser = argparse.ArgumentParser()
    Command(parser)
    return parser.parse_args(argv)


def test_quality_preset_accepted_by_cli():
    args = _parse_qual_score_args(["-i", "f.fq", "-o", "/tmp", "--quality-preset", "novaseq"])
    assert args.quality_preset == "novaseq"


def test_unknown_preset_rejected_by_cli():
    with pytest.raises(SystemExit):
        _parse_qual_score_args(["-i", "f.fq", "-o", "/tmp", "--quality-preset", "bad_preset"])


def test_all_preset_names_accepted_by_cli():
    for name in QUALITY_PRESETS:
        args = _parse_qual_score_args(["-i", "f.fq", "-o", "/tmp", "--quality-preset", name])
        assert args.quality_preset == name