# Integration suite

Runs the simulator for real and checks that its output describes itself.

```bash
NEAT_INTEGRATION=1 pytest tests/integration
```

Without that variable the tests skip, so `pytest tests/` is unaffected. The matrix takes a
couple of minutes; the unit suite takes a second and a half.

## Why this exists

The unit suite pins behaviour that is already understood — a known bug, a known branch. It is
good at that, and it is not what catches the defects this project actually ships. Every one of
the following passed a green suite and was found by generating reads and asking whether the
output was self-consistent:

| defect | what was wrong | suite at the time |
|---|---|---|
| #324 | a `./.` genotype aborted the run | 858 green |
| #325 | every input variant's metadata silently dropped | 858 green |
| #326 | 77 of 26,666 reads placed at the wrong coordinate | 883 green |
| #335 | 85 of 16,000 records displaced by one base | 918 green |
| #337 | one malformed read truncated a 12,000-record BAM at 7,591 | 919 green |

The common shape: a small fraction of records, wrong in a way that exits 0 and produces a
complete, individually well-formed output set. Nothing short of running the thing and checking a
property over every record will find those.

## What is checked

Each invariant lives in `invariants.py` and returns a list of failures rather than asserting, so
one run reports everything wrong with it.

- **Reads are placed where their alignment says.** Walk each record's own POS and CIGAR against
  the reference; a record is misplaced only if it matches substantially better at a nearby
  offset. Sequencing error alone can push the mismatch rate high, so the test is relative rather
  than absolute.
- **CIGARs account for every base.** `M + I + S` must equal the read length.
- **CIGARs are well formed.** No leading or trailing `D`.
- **The golden BAM is sorted and indexable.** Per contig, and the index is actually built —
  the runner indexes it, so an unsorted record aborts the run after the FASTQs are on disk.
- **FASTQ and BAM agree.** Same record count, one quality character per base.
- **The golden VCF is valid.** `bcftools` accepts it, and every INFO key and FILTER value a
  record uses is declared in the header NEAT wrote.
- **Adapter clips are on the right side.** Start of the CIGAR for a reverse read, end for a
  forward one.

## The matrix

The axes are the ones that have actually surfaced defects, not every option NEAT has:

- `read_len` 75 / 150 / 250 — only 250 exposed #335
- `fragment_mean` relative to `read_len` — short inserts make the mates cover the same window and
  remove the headroom a deletion needs
- adapters on/off, and the presets
- paired vs single-ended — single-ended produces no reverse reads at all, so it exercises none of
  the reverse-strand coordinate handling
- `include_vcf` — input variants reach the golden VCF by a different path, and an MNP becomes an
  `UnknownVariant`
- reference composition — uniform, GC-rich, repetitive, multi-contig

`test_input_vcf_and_seeds.py` adds two things the matrix cannot cover: the input-VCF path end to
end, and the same invariants across several seeds. That second one matters because an invariant
holding for one seed may be holding by cancellation rather than by construction — #337 is exactly
that.

## Known failures

`KNOWN_MISPLACED` and `KNOWN_TRAILING_DELETION` in `test_output_invariants.py` record
configurations with an open defect, so the suite is green when nothing *new* is wrong. A nightly
that is permanently red gets ignored. Both entries name the issue; delete one when its issue
closes and the test simply passes.

## Adding a case

Append to `MATRIX`. If a new configuration exposes something, that is the point — file it, and
add it to the known-failures map with the issue number rather than weakening the check.
