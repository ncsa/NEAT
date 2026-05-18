# NEAT v4.4.3
Major performance and memory overhaul focused on making NEAT viable for large
genomes on supercomputing hardware. No user-visible API changes (other than
chunk size now auto-tuning by default).

**Benchmark (ecoli 10× coverage, 4 threads, identical configs):**

| Metric                  | v4.4.2     | v4.4.3     | Improvement     |
|-------------------------|------------|------------|-----------------|
| ecoli SE wall time      | 14:55      | 1:57       | 7.6× faster     |
| ecoli PE wall time      | 14:46      | 1:44       | 8.5× faster     |
| ecoli SE total CPU      | 3,227 s    | 356 s      | 9.1× less       |
| ecoli PE total CPU      | 3,168 s    | 368 s      | 8.6× less       |
| Peak resident memory    | 549 MB     | 175 MB     | 3.1× less       |
| Peak heap (memray)      | 1.27 GB    | 0.32 GB    | 4× less         |
| Per-worker memory       | O(N×cov)   | O(1)       | bounded         |
| `pysam.sort` calls      | 2          | 0          | gone            |
| BAM correctness         | 0.06% dups | strict     | fixed           |

**Versus NEAT 2.1 (single-threaded baseline):**
- SE: 12:28 → 1:57 (6.4× faster, 53% less CPU)
- PE: 20:12 → 1:44 (11.6× faster, 70% less CPU)

**What changed in the hot path:**
- Vectorized error sampling in `get_sequencing_errors` — replaced a ~150-iteration
  per-read Python loop with batched numpy. Eliminated 28M `np.prod` calls per
  185k-read run.
- Vectorized `get_quality_scores` — replaced per-base scalar `rng.normal` with
  one batched call.
- Replaced per-read `PairwiseAligner.align()` in `make_cigar` with a direct
  walker that builds the CIGAR from known error/mutation positions in O(L).
  99% of reads now skip alignment entirely.

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
  threads.

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
