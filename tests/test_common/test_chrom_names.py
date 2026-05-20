"""
Tests for neat/common/chrom_names.py — chromosome-name normalization helpers.
"""
from pathlib import Path

from neat.common.chrom_names import (
    apply_aliases,
    find_aliases,
    load_chrom_aliases,
    prefix_flip_candidates,
)


# ===========================================================================
# prefix_flip_candidates
# ===========================================================================

def test_prefix_flip_strips_chr_prefix():
    assert prefix_flip_candidates("chr1") == {"1"}


def test_prefix_flip_adds_chr_prefix():
    assert prefix_flip_candidates("1") == {"chr1"}


def test_prefix_flip_handles_sex_chroms():
    assert prefix_flip_candidates("chrX") == {"X"}
    assert prefix_flip_candidates("Y") == {"chrY"}


def test_prefix_flip_mitochondrial_chrM_yields_all_three_alternatives():
    """`chrM` should yield {M, MT, chrMT} — all known mt-name conventions."""
    assert prefix_flip_candidates("chrM") == {"M", "MT", "chrMT"}


def test_prefix_flip_mitochondrial_MT_yields_all_three_alternatives():
    assert prefix_flip_candidates("MT") == {"M", "chrM", "chrMT"}


def test_prefix_flip_excludes_self():
    """The input name is never present in its own candidate set."""
    for name in ("chr1", "1", "chrM", "MT", "scaffold_42"):
        assert name not in prefix_flip_candidates(name)


def test_prefix_flip_custom_scaffold_name_only_does_prefix():
    """Non-standard names get the chr-prefix heuristic but no mitochondrial mapping."""
    assert prefix_flip_candidates("scaffold_42") == {"chrscaffold_42"}


# ===========================================================================
# find_aliases
# ===========================================================================

def test_find_aliases_native_match_omitted():
    """Names that already appear in target are skipped."""
    assert find_aliases({"chr1"}, {"chr1", "chr2"}) == {}


def test_find_aliases_prefix_flip_proposed():
    assert find_aliases({"1", "2"}, {"chr1", "chr2"}) == {"1": "chr1", "2": "chr2"}


def test_find_aliases_mitochondrial_proposed():
    assert find_aliases({"MT"}, {"chrM"}) == {"MT": "chrM"}


def test_find_aliases_unmatchable_omitted():
    """A name with no candidate in target is not in the result — caller knows the
    BED is genuinely non-overlapping rather than just misnamed."""
    assert find_aliases({"weird_contig"}, {"chr1"}) == {}


def test_find_aliases_mixed_some_aliased_some_native():
    """A BED can have both naming conventions; aliases only cover mismatches."""
    result = find_aliases({"chr1", "2", "MT"}, {"chr1", "chr2", "chrM"})
    assert result == {"2": "chr2", "MT": "chrM"}


def test_find_aliases_empty_inputs():
    assert find_aliases(set(), {"chr1"}) == {}
    assert find_aliases({"chr1"}, set()) == {}


# ===========================================================================
# load_chrom_aliases
# ===========================================================================

def test_load_chrom_aliases_returns_empty_for_none():
    assert load_chrom_aliases(None) == {}


def test_load_chrom_aliases_parses_tab_separated(tmp_path):
    f = tmp_path / "aliases.tsv"
    f.write_text("1\tchr1\n2\tchr2\nMT\tchrM\n")
    assert load_chrom_aliases(f) == {"1": "chr1", "2": "chr2", "MT": "chrM"}


def test_load_chrom_aliases_skips_comments_and_blanks(tmp_path):
    f = tmp_path / "aliases.tsv"
    f.write_text(
        "# This is a comment\n"
        "\n"
        "1\tchr1\n"
        "# another comment\n"
        "2\tchr2\n"
    )
    assert load_chrom_aliases(f) == {"1": "chr1", "2": "chr2"}


def test_load_chrom_aliases_tolerates_whitespace_separator(tmp_path):
    """Human-typed files often use spaces instead of literal tabs."""
    f = tmp_path / "aliases.tsv"
    f.write_text("1    chr1\n2  chr2\n")
    assert load_chrom_aliases(f) == {"1": "chr1", "2": "chr2"}


def test_load_chrom_aliases_skips_malformed_lines(tmp_path):
    f = tmp_path / "aliases.tsv"
    f.write_text("just_one_token\n1\tchr1\n")
    assert load_chrom_aliases(f) == {"1": "chr1"}


def test_load_chrom_aliases_later_entries_override_earlier(tmp_path):
    """If a source name is listed twice, the last entry wins (documents behavior)."""
    f = tmp_path / "aliases.tsv"
    f.write_text("1\tchr1\n1\tchrONE\n")
    assert load_chrom_aliases(f) == {"1": "chrONE"}


# ===========================================================================
# apply_aliases
# ===========================================================================

def test_apply_aliases_returns_canonical_when_mapped():
    assert apply_aliases("1", {"1": "chr1"}) == "chr1"


def test_apply_aliases_returns_input_when_unmapped():
    assert apply_aliases("chr1", {"1": "chr1"}) == "chr1"


def test_apply_aliases_empty_map_passthrough():
    assert apply_aliases("chr1", {}) == "chr1"
