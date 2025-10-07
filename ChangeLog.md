# NEAT has a new home
NEAT is now a part of the NCSA github and active development will continue here. Please direct issues, comments, and requests to the NCSA issue tracker. Submit pull requests here insead of the old repo.

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

