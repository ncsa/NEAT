# The NEAT Project v4.0
Welcome to the NEAT project, the NExt-generation sequencing Analysis Toolkit, version 4.0. This is our first (beta) release of the newest version of NEAT. There is still lots of work to be done. See the [ChangeLog](ChangeLog.md) for notes.

We are also working on redeveloping NEAT in Rust, a memory and thread safe language that will lend itself well to the way NEAT works, check that out here: https://github.com/ncsa/rusty-neat

Stay tuned over the coming weeks for exciting updates to NEAT, and learn how to [contribute](CONTRIBUTING.md) yourself. If you'd like to use some of our code, no problem! Just review the [license](LICENSE.md), first.

NEAT's read-simulator is a fine-grained read simulator. It simulates real-looking data using models learned from specific datasets. There are several supporting utilities for generating models used for simulation and for comparing the outputs of alignment and variant calling to the golden BAM and golden VCF produced by NEAT.

This is release v4.0 of the software. While it has been tested, it does represent a shift in the software with the introduction of a configuration file. For a stable release using the old command line interface, please see: [NEAT 3.0](https://github.com/ncsa/NEAT/releases/tag/3.3) (or check out older tagged releases)

To cite this work, please use:

> Stephens, Zachary D., Matthew E. Hudson, Liudmila S. Mainzer, Morgan Taschuk, Matthew R. Weber, and Ravishankar K. Iyer. "Simulating next-generation sequencing datasets from empirical mutation and sequencing models." PloS one 11, no. 11 (2016): e0167047.


Table of Contents
=================

  * [neat-genreads](#neat-genreads)
  * [Table of Contents](#table-of-contents)
    * [Requirements](#requirements)
    * [Installation] (#installation)
    * [Usage](#usage)
    * [Functionality](#functionality)
    * [Examples](#examples)
      * [Whole genome simulation](#whole-genome-simulation)
      * [Targeted region simulation](#targeted-region-simulation)
      * [Insert specific variants](#insert-specific-variants)
      * [Single end reads](#single-end-reads)
      * [Large single end reads](#large-single-end-reads)
      * [Parallelizing simulation](#parallelizing-simulation)
  * [Utilities](#utilities)
    * [compute_gc_bias](#computegcbias)
    * [model_fragment_lengths](#modelfraglen)
    * [gen_mut_model](#genmutmodel)
    * [model_sequencing_error](#modelseqerror)
      * [Note on Sensitive Patient Data](#note-on-sensitive-patient-data)


## Requirements (the most up-to-date requirements are found in the environment.yml file)

* Python == 3.10.*
* poetry
* biopython == 1.79
* pkginfo
* matplotlib
* numpy
* pyyaml
* scipy
* pysam
* frozendict

## Installation
To install NEAT, you must create a virtual environment using a tool such as conda. Once activated, you can
use the poetry module in build a wheel file, which can then be pip installed. You will need to run these
commands from within the NEAT directory.

```
> conda env create -f environment.yml -n neat
> conda activate neat
> poetry build
> pip install dist/neat*whl
```

Alternatively, if you wish to work with NEAT in the development environment, you can use poetry install within
the NEAT repo, after creating the conda environment:
```
> conda env create -f environment.yml -n neat
> conda activate neat
> poetry install
```

Test your install by running:
```
> neat --help
```

## Usage
NEAT's core functionality is invoked using the read-simulator command. Here's the simplest invocation of read-simulator using default parameters. This command produces a single ended fastq file with reads of length 101, ploidy 2, coverage 10X, using the default sequencing substitution, GC% bias, and mutation rate models.

Contents of neat_config.yml
```
reference: /path/to/my/genome.fa
```

```
neat read-simulator -c neat_config.yml -o simulated_data
```

The output prefix should not specify a file extension (i.e., .fasta, .fastq, etc.),
as these will be appended by NEAT.

A config file is required. The config is a yml file specifying the input parameters. The following is a brief
description of the potential inputs in the config file. See NEAT/config_template/template_neat_config.yml for a
template config file to copy and use for your runs.

reference: full path to a fasta file to generate reads from
read_len: The length of the reads for the fastq (if using). Integer value, default 101.
coverage: desired coverage value. Float or int, default = 10
ploidy: Desired value for ploidy (# of copies of each chromosome in the organism). Default is 2
paired_ended: If paired-ended reads are desired, set this to True. Setting this to true requires
    either entering values for fragment_mean and fragment_st_dev or entering the path to a
    valid fragment_model.
fragment_mean: Use with paired-ended reads, set a fragment length mean manually
fragment_st_dev: use with paired-ended reads, set the standard deviation of the fragment length dataset

The following values can be set to true or omitted to use defaults. If True, NEAT will produce the file type.
The default is given:

produce_bam: False
produce_vcf: False
produce_fasta: False
produce_fastq: True

error_model: full path to an error model generated by NEAT. Leave empty to use default model
    (default model based on human, sequenced by Illumina)
mutation_model: full path to a mutation model generated by NEAT. Leave empty to use a default
    model (default model based on human data sequenced by Illumina)
fragment_model: full path to fragment length model generate by NEAT. Leave empty to use default model
    (default model based on human data sequenced by Illumina)
gc_model: Full path to model for correlating GC concentration and coverage, produced by NEAT.
    (default model is based on human data, sequenced by Illumina)

partition_mode: by chromosome ("chrom"), or subdivide the chromosomes ("subdivision").
    Note: this feature is not yet fully implemented
threads: The number of threads for NEAT to use.
    Note: this feature is not yet fully implemented
avg_seq_error: average sequencing error rate for the sequencing machine. Use to increase or
    decrease the rate of errors in the reads. Float betwoon 0 and 0.3. Default is set by the error model.
rescale_qualities: rescale the quality scores to reflect the avg_seq_error rate above. Set True to activate if you 
    notice issues with the sequencing error rates in your datatset.
include_vcf: full path to list of variants in vcf format to include in the simulation. These will be inserted as they 
    appear in the input VCF into the final VCF, and the corresponding fastq and bam files, if requested.
target_bed: full path to list of regions in bed format to target. 
    All areas outside these regions will have coverage of 0.
off_target_scalar: manually set the off-target-scalar when using a target bed (if you want to have some percentage of 
    reads from outside the targeted regions. Default is 0. (i.e., setting this to 0.02 would mean off-target areas will 
    have a coverage of ~2% of the total coverage). This is an experimental feature.
discard_bed: full path to a list of regions to discard, in BED format.
mutation_rate: Desired rate of mutation for the dataset. Float between 0.0 and 0.3
    (default is determined by the mutation model)
mutation_bed: full path to a list of regions with a column describing the mutation rate of that region,
    as a float with values between 0 and 0.3. The mutation rate must be in the third column as, e.g., mut_rate=0.00.
no_coverage_bias: Set to true to produce a dataset free of coverage bias
rng_seed: Manually enter a seed for the random number generator. Used for repeating runs. Must be an integer.
min_mutations: Set the minimum number of mutations that NEAT should add, per contig. Default is 0. We recommend setting 
    this to at least one for small chromosomes, so NEAT will produce at least one mutation per contig.
fasta_per_ploid: Produce one fasta per ploid. Default behavior is to produce
    a single fasta showing all variants.                                                                                                                                                                        |

The command line options for NEAT are as follows:

Universal options can be applied to any subfunction. The commands should come before the function name (e.g., neat --log-level DEBUG read-simulator ...), excetp -h or --help, which can appear anywhere in the command.
| Universal Options   | Description                          |
|---------------------|--------------------------------------|
| -h, --help          | Displays usage information           |
| --no-log            | Turn off log file creation           |
| --log-dir LOG_DIR   | Sets the log directory to custom path (default is current working directory |
| --log-name LOG_NAME | Custom name for log file (default is timestamped) |
| --log-level VALUE   | VALUE must be one of [DEBUG, INFO, WARN, WARNING, ERROR] - sets level of log to display |
| --log-detal VALUE   | VALUE must be one of [LOW, MEDIUM, HIGH] - how much info to write for each log record |
| --silent-mode       | Writes logs, but suppresses stdout messages |

read-simulator command line options
| Option              | Description                         |
|---------------------|-------------------------------------|
| -c VALUE, --config VALUE | The VALUE should be the name of the config file to use for this run |
| -o OUTPUT, --output OUTPUT | The path, including filename prefix, to use to write the output files |

## Functionality

![Diagram describing the way that genReads simulates datasets](docs/NEATNEAT.png "Diagram describing the way that gen_reads simulates datasets")

NEAT produces simulated sequencing datasets. It creates FASTQ files with reads sampled from a provided reference genome, using sequencing error rates and mutation rates learned from real sequencing data. The strength of NEAT lies in the ability for the user to customize many sequencing parameters, produce 'golden,' true positive datasets. We are working on expanding the functionality even further to model more species, generate larger variants, model tumor/normal data and more!

Features:

- Simulate single-end and paired-end reads 
- Custom read length
- Can introduce random mutations and/or mutations from a VCF file
  - Supported mutation types include SNPs, indels (of any length), inversions, translocations, duplications
  - Can emulate multi-ploid heterozygosity for SNPs and small indels
- Can simulate targeted sequencing via BED input specifying regions to sample from
- Can accurately simulate large, single-end reads with high indel error rates (PacBio-like) given a model
- Specify simple fragment length model with mean and standard deviation or an empirically learned fragment distribution using utilities/computeFraglen.py
- Simulates quality scores using either the default model or empirically learned quality scores using utilities/fastq_to_qscoreModel.py
- Introduces sequencing substitution errors using either the default model or empirically learned from utilities/
- Accounts for GC% coverage bias using model learned from utilities/computeGC.py
- Output a VCF file with the 'golden' set of true positive variants. These can be compared to bioinformatics workflow output (includes coverage and allele balance information)
- Output a BAM file with the 'golden' set of aligned reads. These indicate where each read originated and how it should be aligned with the reference
- Create paired tumour/normal datasets using characteristics learned from real tumour data
- Parallelization. COMING SOON!

## Examples

The following commands are examples for common types of data to be generated. The simulation uses a reference genome in fasta format to generate reads of 126 bases with default 10X coverage. Outputs paired fastq files, a BAM file and a VCF file. The random variants inserted into the sequence will be present in the VCF and all of the reads will show their proper alignment in the BAM. Unless specified, the simulator will also insert some "sequencing error" -- random variants in some reads that represents false positive results from sequencing.

### Whole genome simulation
Simulate whole genome dataset with random variants inserted according to the default model. 

```
[contents of neat_config.yml]
reference: hg19.fa
read_len: 126
produce_bam: True
produce_vcf: True
paired_ended: True
fragment_mean: 300
fragment_st_dev: 30

neat read-simulator                  \
        -c neat_config.yml          \
        -o /home/me/simulated_reads
```

### Targeted region simulation
Simulate a targeted region of a genome, e.g. exome, with a targeted bed:

```
[contents of neat_config.yml]
reference: hg19.fa
read_len: 126
produce_bam: True
produce_vcf: True
paired_ended: True
fragment_mean: 300
fragment_st_dev: 30
targed_bed: hg19_exome.bed

neat read-simulator                 \
        -c neat_config              \
        -o /home/me/simulated_reads

```

### Insert specific variants
Simulate a whole genome dataset with only the variants in the provided VCF file using -v and setting mutation rate to 0 with -M.

```
[contents of neat_config.yml]
reference: hg19.fa
read_len: 126
produce_bam: True
produce_vcf: True
paired_ended: True
fragment_mean: 300
fragment_st_dev: 30
input_variants: NA12878.vcf
mutation_rate: 0

neat read-simulator                 \
        -c neat_config.yml          \
        -o /home/me/simulated_reads
```

### Single end reads
Simulate single end reads by omitting paired ended options.

```
[contents of neat_config.yml]
reference: hg18.fa
read_len: 126
produce_bam: True
produce_vcf: True

neat read-simulator                 \
        -c neat_config.yml          \
        -o /home/me/simulated_reads
```

### Large single end reads
Simulate PacBio-like reads by providing an error model.
(Not yet implemented in NEAT 4.0)

```
[contents of neat-config.yml]
reference: hg19.fa
read_len: 5000
error_model: errorModel_pacbio.pickle.gz
avg_seq_error: 0.1

neat read-simulator                 \
        -c neat_config.yml          \
        -o /home/me/simulated_reads
```

# Utilities	
Several scripts are distributed with gen_reads that are used to generate the models used for simulation.

## neat compute_gc_bias

Computes GC% coverage bias distribution from sample (bedrolls genomecov) data.
Takes .genomecov files produced by BEDtools genomeCov (with -d option).
(Not yet implemented in NEAT 4.0)

```
bedtools genomecov
        -d                          \
        -ibam normal.bam            \
        -g reference.fa
```

```
neat compute_gc_bias                \
        -r reference.fa             \
        -i genomecovfile            \
        -w [sliding window length]  \
        -o /path/to/prefix
```

## neat model-fraglen

Computes empirical fragment length distribution from sample data.

    neat model-fraglen   \
        -i input.bam            \
        -o /path/to/prefix

and creates fraglen.pickle.gz model in working directory.

## neat gen-mut-model

Takes references genome and VCF file to generate mutation models:

```
neat gen-mut-model reference.fa input_variants.vcf   \
        -o /home/me/models
```

Trinucleotides are identified in the reference genome and the variant file. Frequencies of each trinucleotide transition are calculated and output as a pickle (.p) file.

| Option          | Description                                                                  |
|-----------------|------------------------------------------------------------------------------|
| -o <str>        | Path to output file and prefix. Required.                                    |
| --bed           | Flag that indicates you are using a bed-restricted vcf and fasta (see below) |
| --save-trinuc   | Save trinucleotide counts for reference                                      |
| --human-sample  | Use to skip unnumbered scaffolds in human references                         |
| --skip-common   | Do not save common snps or high mutation areas                               |

## neat model-seq-err

Generates sequencing error model for neat.
This script needs revision, to improve the quality-score model eventually, and to include code to learn sequencing errors from pileup data.

Note that model-seq-err does not allow for sam input, as genSeqErrorModel.py did. If you would like to use
a bam/sam file, use samtools to convert to a fastq, then use the fastq as input.

Note additionally that the file must either be unzipped or bgzipped. If your file is currently gzipped, you can
use samtools' built-in bgzip utility to convert.

gzip -d my_file.fastq.gz
bgzip my_file.fastq

This blocked zip format allows for indexing of the file.

For quality scores, note that using a single number will check quality scores for every number. As this could 
potentially slow down model creation, binned quality scores are advisable.

Soon we will take a samtools mpileup output as input and have some plotting features. 

```
neat model-seq-err                                    \
        -i input_read.fq (.gz)                       \
        -o /path/to/prefix                             \
        -q quality score offset [33]                  \
        -Q maximum quality score [2, 11, 24, 37]      \
        -m maximum number of reads to process [all]   \
```

Please note that -i2 can be used in place of -i to produce paired data.

## neat plot_mutation_model

Performs plotting and comparison of mutation models generated from genMutModel.py (Not yet implemented in NEAT 4.0).

```
neat plot_mutation_model                                                \
        -i model1.pickle.gz [model2.pickle.gz] [model3.pickle.gz]...    \
        -l legend_label1 [legend_label2] [legend_label3]...             \
        -o path/to/pdf_plot_prefix
```

## neat vcf_compare

Tool for comparing VCF files (Not yet implemented in NEAT 4.0).

```
neat vcf_compare
        -r <ref.fa>        * Reference Fasta                           \
        -g <golden.vcf>    * Golden VCF                                \
        -w <workflow.vcf>  * Workflow VCF                              \
        -o <prefix>        * Output Prefix                             \
        -m <track.bed>     Mappability Track                           \
        -M <int>           Maptrack Min Len                            \
        -t <regions.bed>   Targetted Regions                           \
        -T <int>           Min Region Len                              \
        -c <int>           Coverage Filter Threshold [15]              \
        -a <float>         Allele Freq Filter Threshold [0.3]          \
        --vcf-out          Output Match/FN/FP variants [False]         \
        --no-plot          No plotting [False]                         \
        --incl-homs        Include homozygous ref calls [False]        \
        --incl-fail        Include calls that failed filters [False]   \
        --fast             No equivalent variant detection [False]
```
Mappability track examples: https://github.com/zstephens/neat-repeat/tree/master/example_mappabilityTracks

### Note on Sensitive Patient Data
ICGC's "Access Controlled Data" documentation can be found at <a href = https://docs.icgc.org/portal/access/ target="_blank">https://docs.icgc.org/portal/access/</a>. To have access to controlled germline data, a DACO must be submitted. Open tier data can be obtained without a DACO, but germline alleles that do not match the reference genome are masked and replaced with the reference allele. Controlled data includes unmasked germline alleles.
