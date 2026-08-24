# NEAT v4.7.0

Optional 3' sequencing-adapter readthrough, ported from NEAT's Rust sister project
`eidolon` (issue #125 there), plus the adapter-free short-insert control arm that
makes it interpretable.

A sequencer runs a fixed number of cycles regardless of insert size, so a fragment
whose insert is shorter than the read length gets read through into the adapter
ligated at its far end. NEAT could not produce this. `fragment_mean` was required to
be at least `read_len`, and both fragment samplers discarded anything below
`read_len + 10` and resampled, which silently left-truncated the insert-size
distribution the user asked for: with `read_len: 150`, `fragment_mean: 300`,
`fragment_st_dev: 100`, about 8% of the distribution was dropped and the surviving
mean was ~16 bp too long. Every read was genomic by construction, so a trimmer run
against NEAT output reported 0% adapter whether or not it was correctly configured —
a QC step that could not fail — and short-insert library types (cfDNA, FFPE and other
degraded DNA, ancient DNA, small RNA, ATAC-seq nucleosome-free fragments) could not
be simulated at all.

New config options:

- `adapters` — model 3' adapter readthrough. Inserts shorter than `read_len` are kept
  rather than resampled, and each read is padded back out to `read_len` at its 3' end
  with adapter sequence (read 1 gets `adapter_r1`, read 2 gets `adapter_r2`, appended
  after read 2 is reverse-complemented so it lands in read orientation). Adapter bases
  take quality scores and substitution errors from the run's models — substitutions
  only, since an indel there would break the exact read length the padding provides —
  are soft-clipped in the golden BAM, and carry no variants.
- `adapter_preset` — `truseq` (default), `nextera`, or `custom`.
- `adapter_r1` / `adapter_r2` — explicit 5'->3' uppercase A/C/G/T sequences, required
  when `adapter_preset` is `custom`.
- `keep_short_fragments` — keep short inserts as plain insert-length genomic reads with
  no adapter padding. Implied by `adapters`; set on its own it is the isolation control
  that separates the adapter effect from the reduced-coverage effect of short inserts.

Enabling either feature relaxes the `fragment_mean >= read_len` check to a warning, so
short-insert libraries can be configured. Inserts below 25 bp are still resampled, which
matches real size selection: a read that is almost entirely adapter carries no usable signal.

`FragmentLengthModel.generate_fragments` no longer splices a fixed set of tiny lengths
(10, 11, 12, 13, 14, 28, 31 bp) into every batch. Those values were anti-infinite-loop
padding rather than draws from anyone's distribution, and the ordinary `read_len` fragment
floor had been hiding them; the lower floor that short-insert runs use admits the 28 and
31 bp entries, which would have put a spike of artificial, adapter-heavy inserts into the
output. The samplers now bound their own retries instead, and report an under-covered
region rather than spinning when a fragment model cannot clear the floor. Fragment lengths
for a given seed differ from previous releases as a result.

A deletion falling inside a short insert is now applied rather than dropped. Deletion
headroom is normally literal — reference just past the read, which a deletion pulls in to
keep the read at `read_len` — but a short insert is the whole molecule and has none. It is
therefore budgeted instead: the deletion shortens the genomic portion, and the adapter tail
grows by the same amount so the read still reaches `read_len` (with `keep_short_fragments`
alone there is nothing to backfill with, so the read is emitted shorter, as a shorter
sequenced molecule should be). Previously the zero headroom made every guard keyed on
padding skip the deletion wherever it fell, so a short-insert read came back identical to
the unmutated reference while the golden VCF still claimed the variant — on a 4000-read
paired run with 300 deletions, 60 reads carried one where 1371 now do. A deletion that
genuinely cannot be applied is now logged at debug level rather than dropped silently.

One known limitation remains. A deletion at the extreme 3' edge of a reverse short-insert
read lands on the CIGAR's leading edge once SEQ is flipped to reference-forward, where the
alignment absorbs it into POS instead of emitting a leading `D`. The read is correct and the
record well formed; the placement of that single edge base is not.

Applying a deletion *error* now keeps the quality score of its anchor base. The anchor survives
in the read -- the alternate allele is that base, VCF-style -- but its score was dropped, leaving
the quality array one shorter than the sequence per deletion error. The defect is old and was
masked twice over: a full-length read draws its quality array across the padded reference segment
and trims to the read length afterwards, which absorbed the loss, and a short insert could never
reach it while zero padding made every deletion be skipped. Fixing short-insert deletions removes
both masks, and the malformed records that result accumulate into a golden BAM that `samtools
index` rejects outright.

Because SAM stores QUAL in the same orientation as SEQ, a reverse read's quality array is
now reversed alongside its sequence when the golden BAM record is written. Flipping only the
sequence paired every base with another base's score — invisible on a full-length read with
a flat quality profile, but on a short-insert read the adapter tail's scores landed on the
genomic prefix.

Because SAM stores SEQ in reference-forward orientation, a reverse read's 3' adapter tail
is written as a *leading* soft clip in the CIGAR and a forward read's as a trailing one.

Note on coverage: short inserts carry fewer genomic bases per read and their mates overlap,
so realized depth falls below the requested `coverage`. This matches `eidolon`'s behavior and
is documented rather than compensated for, so the two simulators stay comparable arm-for-arm.

Both features default to off, and output is unchanged when their keys are omitted —
verified byte-for-byte on FASTQ and on all BAM alignment records for a fixed seed.

# NEAT v4.6.2

Bug fixes for input VCF (`include_vcf`) handling. Both issues were reported against
the v4.6.1 bioconda build.

- A no-call genotype (`./.`) in an input VCF aborted the run with
  `ValueError: invalid literal for int() with base 10: '.'` before any read was
  generated, leaving an empty output directory (#324). A no-call is valid VCF 4.2
  and is what several callers emit for a record that is present but unassigned to a
  genotype. It is now treated the way a record whose FORMAT has no `GT` at all was
  already treated: a genotype is generated with `pick_ploids`, with a warning and a
  tally, and the sample column is rewritten so the reads and the golden VCF agree on
  the genotype used. Other FORMAT subfields are preserved.

- Input variants were constructed with `kwargs=data` against a `**kwargs`
  signature, so every one of them carried `metadata = {'kwargs': {...}}` and every
  metadata lookup missed. This dropped input `ID`/`FILTER`/`INFO`/`FORMAT` from the
  golden VCF silently, and raised `KeyError: 'REF'` from
  `UnknownVariant.get_ref_len` as soon as a read carried a multi-nucleotide variant
  — coverage-dependent, and after the FASTQs had already been written, so the run
  left a partial output set (#325). `UnknownVariant` also never set `self.alt`.

- Restoring that metadata meant it began reaching the output, which exposed two
  latent defects that would otherwise have broken previously-working runs,
  including the H1N1 example in `data/`:

  - The golden VCF carried input `ID`/`FILTER`/`INFO`/`FORMAT` values without the
    matching `##INFO`/`##FILTER`/`##FORMAT` declarations, which is not valid VCF.
    `bcftools sort` fails on the first such record and the run then died with a
    `FileNotFoundError` on its temp file. Declarations from the input header are now
    carried over, and permissive ones are synthesised for keys the input itself left
    undeclared (`data/H1N1.vcf` uses `AF`, `PP` and `DP` while declaring none).
  - Records whose FORMAT was `.` produced `GT:.`, naming a FORMAT key `.`. A missing
    FORMAT is now replaced rather than prepended to.

- `alt_count` split the ALT column on `;` rather than `,`. ALT alleles are
  comma-separated, so the count was always 1 and `pick_ploids` could never draw the
  second alternate of a multiallelic record.

# NEAT v4.6.1

Internal cleanup follow-up to the v4.6.0 N-handling work.

- Removed a dead `has_r2` guard in `generate_reads._filter_n_regions`. In the
  paired branch read2 is always a real, in-bounds window (`(e - read_len, e)`
  with `e >= read_len + 10`), so the guard was always true; the `(0, 0)`
  placeholder mate only occurs in single-ended mode, which never reaches that
  branch. Behavior is unchanged (covered by the paired N-filter tests).

No user-facing behavior changes.

# NEAT v4.6.0

Realistic handling of reference `N` (unknown) bases.

Previously every `N` in a read's reference window was filled with a low-quality
`TTAGGG` human-telomere repeat (`Read.convert_masking`), and reads were never
excluded from `N` regions. This fabricated the single most mismap-prone sequence
in the genome — the telomere repeat occurs at every chromosome end — at default
mapping quality, inside assembly gaps (centromeres, telomeres, heterochromatin)
where real WGS produces essentially zero coverage. It was also biologically wrong
for the many non-human genomes NEAT simulates.

NEAT now models `N` regions the way real data behaves, via two complementary steps
(both on by default):

- **Placement exclusion:** fragments whose read window is at least `n_max_fraction`
  (default 0.5) `N` are dropped at sampling time, so true assembly gaps get ~zero
  coverage. Coverage of callable sequence is unaffected.
- **Literal-N edges:** an `N` that remains at the edge of a surviving read is emitted
  as a literal `N` base call at floor quality, rather than fabricated sequence.
  Aligners and callers treat `N` as no-information, not as a mismatch.

New config options:

- `n_handling` — `exclude` (new default, the behavior above) or `telomere` (the
  legacy `TTAGGG` masking, retained for reproducing older runs).
- `n_max_fraction` — `N`-fraction at or above which a read is dropped under
  `exclude` (range 0.0–1.0, default 0.5).

Note: under the default `exclude` policy this changes the reads produced for any
reference containing `N`. Set `n_handling: telomere` to reproduce pre-4.6.0 output.

# NEAT v4.5.3

Handle IUPAC ambiguity codes in the reference (issue #291).

Reference assemblies such as GRCh38 carry IUPAC ambiguity codes (R, Y, S, W,
K, M, B, D, H, V) alongside A/C/G/T/N. NEAT previously had no handling: such a
base could survive reference loading and then crash a run with a `KeyError`
when it landed at a sequencing-error site (the substitution model only indexes
A/C/G/T).

- Ambiguity codes are now resolved to a concrete base — one of the bases the
  code represents, chosen with the run's seeded RNG — when the reference is
  loaded. Every downstream consumer (reads, BAM, golden VCF, error and
  trinucleotide lookups) therefore only ever sees A/C/G/T/N, and a count of
  resolved bases is logged per contig. `N` keeps its existing low-quality
  masking behavior.
- Defense-in-depth: the sequencing-error model now skips a non-A/C/G/T
  reference base instead of raising, so no unmapped base can crash a run.

# NEAT v4.5.2

Fixed coverage handling at low and fractional depths (issue #242).

Previously no reads were generated at `coverage=1`, and fractional coverage was
rejected outright: `coverage` was typed as an integer in the config schema, and
the per-chunk read count was rounded with `math.ceil`.

- `coverage` is now a float with a positive lower bound, so fractional depths
  (e.g. `0.5` for low-pass simulation) are accepted and zero/negative is
  rejected.
- Each chunk's read count is kept as a float expectation and rounded
  stochastically (`E[reads] == target`) rather than with `ceil`. Rounding every
  chunk up systematically over-covered — negligible at high depth but dominant
  at low/fractional coverage. The run's seeded RNG is used, so output remains
  reproducible. Note: same-seed output now differs from prior versions.
- A chunk that legitimately rounds to zero reads no longer crashes the uniform
  sampler.

# NEAT v4.5.1

Removed the `cleanup_splits` and `reuse_splits` config options.

`reuse_splits` was fully broken (it raised `FileNotFoundError` unconditionally
regardless of whether the splits directory existed). `cleanup_splits` existed
solely to support `reuse_splits`; without it the option has no value — the
splits directory always lives in a `TemporaryDirectory` that is cleaned up
automatically when the run exits.

Both keys are now in `DEPRECATED_KEYS`: configs that still include them
receive a one-line deprecation warning and continue parsing cleanly, so no
user config breaks. The README examples and all internal test fixtures have
been updated to drop these keys.

# NEAT v4.5.0

New `neat compare-vcfs` subcommand: compares a downstream variant caller's VCF
against the NEAT-simulated truth VCF and attributes each false negative to the
simulator's own configuration (mutation bed, target bed, simulated contigs).
Issue #297.

Variant equivalence (multi-allelic normalization, haplotype-level matching) is
delegated to Illumina's `hap.py`. NEAT adds the false-negative attribution
layer: each FN is tagged with one or more of `outside_simulated_contigs`,
`outside_mutation_bed`, `outside_target_bed`, or `unknown`.

**Outputs (in `--output-dir`):**

- `comparison_summary.json` — counts, precision/recall/F1, per-reason FN totals
- `comparison_summary.txt` — human-readable rollup
- `FN_with_reasons.vcf` — hap.py's FN records with an added `NEAT_REASON` INFO tag
- `happy.vcf.gz` and siblings — preserved hap.py output
- `fn_attribution.png` — optional bar chart, written when `--plot` is set

**Prerequisite artifact:** Every `neat read-simulator` run now writes a small
`simulation_summary.json` alongside its other outputs, capturing the run's
config echo (coverage, read length, paired-ended, BED paths, contigs simulated)
and delivered counts (total reads, total variants, per-contig variants). The
`compare-vcfs` wrapper reads this file to drive attribution.

**External dependency:** `hap.py` is required at runtime. NEAT does not bundle
it; install via `conda create -n hap_py_env -c bioconda -c conda-forge hap.py`
and pass the absolute path via `--happy-bin`, or put it on `$PATH`. Without
hap.py, the command exits with an install hint.

**Chromosome-name handling:** `compare-vcfs` detects when a BED's chrom names
don't overlap the reference's (e.g., BED uses `1`/`MT` while reference uses
`chr1`/`chrM`) and writes a warning into `comparison_summary.json` suggesting
an alias mapping. Users can apply the mapping via a new `--chrom-aliases TSV`
flag. NEAT does not auto-normalize — silent renaming would mask real bugs.
`simulation_summary.json` now also records `delivered.reference_contigs` (the
full FASTA contig set) alongside `contigs_simulated`.

**Not in this release** (deferred to follow-up issues): full per-region
simulation telemetry (per-chunk coverage, GC-bias map, error rates by position)
for richer FN attribution, and SV-comparison support. The simulator's silent
"skip BED chroms not in reference" warning will be promoted to a fatal error
in a future release with an opt-in `chrom_aliases` config key.

# NEAT v4.4.4
Follow-up release on top of v4.4.3 bundling three lines of work: another perf
pass over the remaining single-thread hot paths (variant overlap checks,
trinucleotide bias mapping, BAM record encoding), a coverage-semantics fix
for the GC bias model added earlier in the 4.4 line, and removal of the
`parallel_mode` config key. The perf changes are behavior-preserving for a
given seed; the GC bias fix changes output bytes when a non-uniform `gc_model`
is configured (a typical model with mean weight ~0.7 previously delivered
~70 % of the requested coverage; v4.4.4 delivers the full coverage).

**Benchmark (ecoli 10× coverage, 4 threads, identical configs):**

| Metric                       | v4.4.3   | v4.4.4   | Improvement     |
|------------------------------|----------|----------|-----------------|
| ecoli SE wall time           | 1:35     | **1:07** | 1.42× faster    |
| ecoli PE wall time           | 1:34     | **1:06** | 1.43× faster    |
| ecoli SE total CPU           | 331 s    | 244 s    | 26% less        |
| ecoli PE total CPU           | 338 s    | 248 s    | 27% less        |
| Worker CPU utilization (PE)  | 359%/400 | **378%** | 94.5% of theory |
| Peak resident memory         | 175 MB   | 175 MB   | unchanged       |

Cumulative vs the v4.4.2 baseline at the start of this performance arc:
SE 14:55 → 1:07 (**13.3× faster**), PE 14:46 → 1:06 (**13.4× faster**).

**Hot-path fixes:**

- `ContigVariants.check_if_del` / `check_if_ins` — replace O(N) linear scans
  of `all_dels` / `all_ins` with O(1) lookups into per-position bucket dicts.
  Each indel is indexed at every reference position it covers when added;
  total memory cost is O(sum of indel lengths). Single-thread profile on
  14 k calls: `check_if_del` 3.62 s → 0.011 s (325×); `check_if_ins` 3.70 s →
  0.007 s (530×). Removes O(N²) scaling that would have hurt larger genomes.

- `variant_models.map_local_trinuc_bias` — replace 64 regex scans per call
  (one `re.finditer` per trinucleotide) with a single vectorized numpy
  expression. The cache check is also tightened from `==` (O(N) byte
  compare) to `is` (O(1) identity), since callers pass different slice
  objects each time anyway. Single-thread profile on 14 k calls: 13.5 s →
  0.81 s (16×). The entire `generate_variants` call dropped from 24.9 s
  cumulative to 4.2 s.

- `OutputFileWriter.write_bam_record` — vectorize BAM record encoding:
  - Sequence encoding replaces the per-byte Python loop (28 M
    `Seq.__getitem__` calls in the old code) with a 256-entry numpy lookup
    plus a single `((codes[0::2] << 4) | codes[1::2]).tobytes()`.
  - Quality encoding replaces `"".join([chr(n) for n in quality_array])`
    with one `np.asarray(...).astype(np.uint8).tobytes()` call.
  - CIGAR ops are now packed in one `struct.pack(f'<{n}I', ...)` instead of
    one pack-per-op + `bytearray.extend`. The CIGAR string is parsed in a
    single linear pass (the previous code did `re.split` + `re.findall`).
  - The fixed-size BAM record header is now one
    `struct.pack('<iiiIIiiii', ...)` instead of nine separate pack calls
    glued with `+`.
  Single-thread profile: self time 12.0 s → 3.2 s, cumtime 58.1 s → 8.1 s
  (7.2× cumulative). This was the single biggest hot path of the release.

The vectorized sequence and CIGAR encodings produce byte-identical BAM records
to the prior implementation (verified by `pysam.AlignmentFile` read-back and
sort-order checks on the 4-thread ecoli benchmark).

**Config cleanup:** `parallel_mode` has been removed from the YAML schema.
The splitting strategy is now derived from `threads` (single-chunk-per-
contig when `threads == 1`, size-based chunking when `threads > 1`).
Existing configs with `parallel_mode: ...` keep parsing — the key is
silently ignored with a one-line deprecation warning. The option never
actually let users change effective behavior under the old logic (it was
force-overridden to `'contig'` when `threads == 1`, and `'size'` was a
strict superset for `threads > 1`), so this is API surface cleanup rather
than a behavior change.

**GC bias coverage fix.** When a non-uniform `gc_model` was configured,
`cover_dataset` had been scaling per-chunk reads by `chunk_mean / max_weight`.
That made `coverage` the *peak* coverage at the GC bias maximum, undercounting
total reads by `mean_weight_genome / max_weight` (a typical model with
mean ~0.7, max 1.0 delivered ~70 % of the requested coverage). The user-facing
config (`template_neat_config.yml:18`) documents `coverage` as "average coverage
for the entire genome," so this was a documentation/behavior mismatch.

The fix:

- New `compute_genome_wide_gc_mean_weight(reference, gc_model)` helper
  (`neat/model_gc_bias/utils.py`) — single-pass vectorized scan over every
  length-`window_size` window in the reference, reusing the same prefix-sum
  windowing as `cover_dataset`. Uniform models short-circuit without scanning.
- `read_simulator_runner` precomputes this once and stores it on `options` so
  all worker chunks see the same denominator (no redundant per-chunk genome
  scans).
- `cover_dataset` now scales by `chunk_mean / global_mean`, which averages to
  1.0 across chunks and preserves `options.coverage` as the genome-wide average.
- When `cover_dataset` is called directly outside the runner (e.g., unit tests
  with an in-memory `SeqRecord`), `gc_global_mean_weight` is unset and the
  scaling falls back to chunk-local mean — i.e., no rescaling, which is the
  correct single-chunk behavior.

The stale `TODO use gc bias to skew this number. Calculate at the runner level.`
at `generate_reads.py:72` is removed; that work is now complete.

All 603 unit/integration tests pass (up from 599 with three new
`test_compute_genome_wide_gc_mean_*` unit tests plus
`test_gc_bias_preserves_average_coverage`, which runs the full pipeline against
a half-AT / half-GC reference under a non-uniform model and asserts delivered
reads are within 20 % of `coverage × genome_length / read_len`).

# NEAT v4.4.3
Major performance and memory overhaul focused on making NEAT viable for large
genomes on supercomputing hardware. No user-visible API changes (other than
chunk size now auto-tuning by default).

**Benchmark (ecoli 10× coverage, 4 threads, identical configs):**

| Metric                  | v4.4.2     | v4.4.3     | Improvement     |
|-------------------------|------------|------------|-----------------|
| ecoli SE wall time      | 14:55      | 1:35       | 9.4× faster     |
| ecoli PE wall time      | 14:46      | 1:35       | 9.4× faster     |
| ecoli SE total CPU      | 3,227 s    | 331 s      | 9.7× less       |
| ecoli PE total CPU      | 3,168 s    | 338 s      | 9.4× less       |
| Peak resident memory    | 549 MB     | 175 MB     | 3.1× less       |
| Peak heap (memray)      | 1.27 GB    | 0.32 GB    | 4× less         |
| Per-worker memory       | O(N×cov)   | O(1)       | bounded         |
| `pysam.sort` calls      | 2          | 0          | gone            |
| BAM correctness         | 0.06% dups | strict     | fixed           |

**Versus NEAT 2.1 (single-threaded baseline):**
- SE: 12:28 → 1:35 (7.9× faster, 56% less CPU)
- PE: 20:12 → 1:35 (12.8× faster, 72% less CPU)

**Scale-test (c_elegans 10× coverage, 4 threads, 100 Mb genome — ~7× the
ecoli reference):**

| Metric                  | Value                                                   |
|-------------------------|---------------------------------------------------------|
| Wall time               | 19:16                                                   |
| Total CPU               | 4,085 s                                                 |
| Peak resident memory    | 304 MB                                                  |
| BAM records             | 6,685,764                                               |
| BAM sort violations     | 0                                                       |
| Stitch step (parallel)  | 5.3 s                                                   |
| Auto-tuned chunk size   | 3.1 Mb (35 chunks)                                      |

Scaling behavior vs ecoli is ~linear in genome size as expected. The stitch
step is bounded by raw disk I/O via `pysam.cat`, so it stays at single-digit
seconds even as the BAM grows. Per-worker peak RSS is 304 MB ÷ 4 ≈ 76 MB,
which is the reference segment + models — independent of chunk size and
coverage.

**What changed in the hot path:**
- Vectorized error sampling in `get_sequencing_errors` — replaced a ~150-iteration
  per-read Python loop with batched numpy. Eliminated 28M `np.prod` calls per
  185k-read run.
- Vectorized `get_quality_scores` — replaced per-base scalar `rng.normal` with
  one batched call.
- Replaced per-read `PairwiseAligner.align()` in `make_cigar` with a direct
  walker that builds the CIGAR from known error/mutation positions in O(L).
  99% of reads now skip alignment entirely.
- Rewrote `apply_errors` as a single ascending-position pass. The previous
  implementation did one `np.concatenate` and one `MutableSeq` slice/concat
  per error — quadratic in errors-per-read. The new pass is linear regardless
  of error count.
- Removed redundant `deepcopy(self.reference_segment)` calls in
  `convert_masking` and `finalize_read_and_write`. Biopython `Seq` is
  immutable; the downstream operations make their own working copies.

**What changed in the I/O path:**
- Removed both `pysam.sort` calls. Per-worker BAMs are emitted coordinate-sorted
  by construction; `pysam.merge` of sorted inputs already produces sorted output.
  The final sort allocated a 1 GB buffer that dominated peak memory.
- Replaced `pysam.merge` with `pysam.cat` for the final stitch. cat does a raw
  BGZF concatenation (no decompression / re-encode), bounded by raw disk I/O
  instead of BGZF rate. At human-30× scale this is the difference between a
  multi-hour stitch and a multi-minute one.
- Each chunk now owns a non-overlapping reference range for read1 placement
  (`responsibility_length`), enabling the cat-based stitch and eliminating
  ~0.06% over-coverage in chunk-overlap regions.
- Streamed FASTQ and BAM records directly to output during read generation.
  Workers no longer accumulate `reads_to_write` — per-worker memory is now
  bounded by reference segment + models, not by chunk size × coverage.
- Stitch steps (FASTQ concat, VCF dedup, BAM cat) now run concurrently in
  threads. On a single-disk system the wall is bounded by the BAM cat alone;
  on parallel filesystems the overlap is more pronounced.
- FASTQ stitch is now byte-level: per-chunk gzip streams are concatenated
  without decompression / re-encode (concatenated gzip streams form a valid
  gzip file per the spec).

**Defaults and ergonomics:**
- `parallel_block_size` now auto-tunes from genome length and thread count
  (target: ~8 chunks per thread). For small bacterial genomes this matches the
  old hardcoded 500 kb; for human-scale genomes it produces ~6 Mb chunks
  instead of ~500 kb, dramatically reducing stitch overhead. Specify the option
  explicitly to override.
- FASTQ output is no longer shuffled; reads come out in the natural sampling
  order. Pipe through `seqkit shuffle` if you need a uniform shuffle (documented
  in README).
- Added a "Multi-node deployment on HPC clusters" section to the README
  showing a SLURM array-job pattern for whole-genome simulation across nodes.

**Caveats:**
- Several of the vectorization fixes change how the PRNG stream is consumed.
  Same seed will produce statistically equivalent reads, but not bit-identical
  to v4.4.2. Re-baseline any regression tests that compared exact output.

# NEAT v4.4.2
- Added GC bias modeling to generate reads and a function to create a GC bias model from real data.
- Added improvements and efficiency upgrades to generate-reads.
- Fixed version number bug.

# NEAT v4.4.1
- Added `readme = "README.md"` to `pyproject.toml` so the project description appears correctly on PyPI.

# NEAT v4.4
- Official release of NEAT 4.0. Represents major contributions from NCSA and beyond.
- Added parallel processing support, making NEAT production-ready for large genomes.
- Includes more options for quality score modeling, a bacterial genome wrapper, additional tests, and performance improvements.
- Fixed inverted indel gate condition in the sequencing error model: indel errors were never generated due to a `>` vs `<=` logic error.
- Fixed deletion blacklist bug: the anchor position of a deletion was incorrectly blacklisted, causing the deletion to remove itself during error cleanup.
- Fixed quality array float promotion: applying a deletion error produced an empty quality slice (`np.array([])`), which defaulted to `float64` and broke quality string encoding.
- Fixed fallback infinite loop in error selection when all quality scores are uniform.
- Fixed off-by-one in the main error selection loop.
- Added `errors_per_read` pre-calculation in the runner for more accurate per-contig error budgeting proportional to coverage and contig length.
- Fixed trinucleotide context slice off-by-one in variant generation.
- Removed debug `import pdb` statements from production code.
- Replaced debug print sentinels in the mutation model with proper log warnings.
- Expanded test coverage for the error model, runner, and read simulator.

# NEAT v4.3.6
- Multiple bug fixes, fixes to outputs. See release for full notes.

# NEAT v4.3.5
- An improvement rather than a bug fix this time. We moved vcf processing into the threaded portion, as our speeds were better than single threaded, but very slow on the vcf writing portion. This sped things up considerably, so we tested and confirmed that it is working as desired and are updating to a new version with improved VCF production in parallel mode.

# NEAT 4.3.4
- More fixes to bam creation, finalizing parallel code

# NEAT 4.3.3
- Bug fixes, but inadvertently pushed a non-working branch. Will remove.

# NEAT v4.3.2
- Bug fixes for parallel processing, which was causing some of the headers to be printed incorrectly. To fix that, we had to rewrite a bunch of the code and integrate parallelism more directly into NEAT.

# NEAT v4.3.1
- Bug fixes (see issue #160) having to do with output files.

# NEAT v4.3.1
- Updated parallel module to integrate it into the code more fluidly. We also updated the options section to revise the process and allow for copying of options objects for parallelism run.

# NEAT v4.3
- Added a parallelization module to run NEAT in parallel. We expect this to speed up times. Please let us know if it works for you!

# NEAT v4.2.X
- Several bug fixes
- Removed bin scoring (needs to be updated, and we didn't have time, but we'll bring it back)
- Incremental improvements

# NEAT v4.2
- After several bug fixes that constituted release 4.1 and some minor releases, we are ready to release an overhauled vesion of NEAT 4.0.
- Removed GC bias - it had little to no effect and made implementation nearly impossible
- Removed fasta creation - we had tweaked this a bit but never got any feedback. It may come back if requested.
- Improvements/fixes/full implementations of:
  - heterozygosity
  - read creation (now with more reads!)
  - bam alignment/creation
  - bed tool incorporation

-Updated "master" branch to "main." - please update your repo accordingly
# NEAT v4.0
- Rewritten the models. Models generated on old versions of NEAT will have to be redone, due to the restructuring of the codebase. These new models should be smaller and more efficient. We have replicated the previous default models in the new style. There is no straightforward way to convert between these, unfortuantely.
- More flexibility. NEAT should be able to handle a wider variety of ploidies, read lengths, and chromosome sizes.
- Binned quality scores. NEAT now assumes binned quality scores.
- File output flexibility. NEAT now can generate each file type independently. Producing a VCF no longer requires NEAT to save all the caluclations of formulating a BAM, which should make the process more efficient.
- Speedups and improvements
- Across the board improvements to code documentation and style

# NEAT v4.0 - beta
- Rewritten codebase to be modular and pip installable, using the poetry package
- Working on optimizing code and implementing final features

# NEAT v3.1

Bug Fixes:
- Fixed an issue where NEAT could not read in gzipped vcf files
- Fixed a bug where NEAT was assuming a ploidy of 2
- Fixed a bug that was causing NEAT to produce exactly the same reverse strands as forward strands in paired ended mode
- Fixed a bug when reading in some mutation models.
- Added gzipping to all mutation models to maximize their compression
- Fixed some issues with vcf_compare

# NEAT v3.0
- NEAT gen_reads now runs in Python 3 exclusively. The previous, Python 2 version is stored in the repo as v2.0, but will not be undergoing active development.
- Converted sequence objects to Biopython Sequence objects to take advantage of the Biopython library
- Converted cigar strings to lists. Now there are simply a functions that convert a generic cigar string to a list and vice versa.
- Tried to take better advantage of some Biopython libraries, such as for parsing fastas.
- For now, we've eliminated the "Jobs" option and the merge jobs function. We plan to implement multi-threading instead as a way to speed up NEAT's simulation process.
- Added a basic bacterial wrapper that will simulate multiple generations of bacteria based on an input fasta, mutate them and then produce the fastqs, bams, and vcfs for the resultant bacterial population.

For improvements, we have a lot of great ideas for general improvements aimed at better simulating bacteria, but we believe this same improvements will have applications in other species as well. 
- Multiploidy - all right this has nothing to do with bacteria specifically, but it is a feature we would like to implement into gen_reads.
- Structural Variants - model large scale structural variants with an eye toward intergenic SVs.
- Transposable Elements - model transposons within the sequence
- Repeat regions - will bring a variety of interesting applications

# NEAT has a new home
NEAT is now a part of the NCSA github and active development will continue here. Please direct issues, comments, and requests to the NCSA issue tracker. Submit pull requests here insead of the old repo.
