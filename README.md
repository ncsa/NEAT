# The NEAT Project v3.2
Welcome to the NEAT project, the NExt-generation sequencing Analysis Toolkit, version 3.2. Neat has now been updated with Python 3, and is moving toward PEP8 standards. There is still lots of work to be done. See the [ChangeLog](ChangeLog.md) for notes.

Stay tuned over the coming weeks for exciting updates to NEAT, and learn how to [contribute](CONTRIBUTING.md) yourself. If you'd like to use some of our code, no problem! Just review the [license](LICENSE.md), first.

NEAT-gen_reads is a fine-grained read simulator. GenReads simulates real-looking data using models learned from specific datasets. There are several supporting utilities for generating models used for simulation.

This is an in-progress v3.2 of the software. For a stable release of the previous repo, please see: [genReads1](https://github.com/zstephens/genReads1) (or check out our v2.0 tag)

To cite this work, please use:

> Stephens, Zachary D., Matthew E. Hudson, Liudmila S. Mainzer, Morgan Taschuk, Matthew R. Weber, and Ravishankar K. Iyer. "Simulating next-generation sequencing datasets from empirical mutation and sequencing models." PloS one 11, no. 11 (2016): e0167047.


Table of Contents
=================

  * [neat-genreads](#neat-genreads)
  * [Table of Contents](#table-of-contents)
    * [Requirements](#requirements)
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
    * [compute_gc.py](#computegcpy)
    * [compute_fraglen.py](#computefraglenpy)
    * [generate_mutation_model.py](#genmutmodelpy)

    * [genSeqErrorModel.py](#genseqerrormodelpy)
    * [plot_mut_model.py](#plotmutmodelpy)
    * [vcf_compare_OLD.py](#vcf_compare_oldpy)
      * [Note on Sensitive Patient Data](#note-on-sensitive-patient-data)


## Requirements

* Python >= 3.6
* biopython == 1.79
* matplotlib >= 3.3.4 (optional, for plotting utilities)
* matplotlib_venn >= 0.11.6 (optional, for plotting utilities)
* pandas >= 1.2.1
* numpy >= 1.22.2
* pysam >= 0.16.0.1

## Notes
Some utilities require outside tools. For example, generate_mutation_model.py requires the user have access to BedTools
if you want to restrict the VCF by a bed file. See specific tools for more information.
## Usage
Here's the simplest invocation of genReads using default parameters. This command produces a single ended fastq file with reads of length 101, ploidy 2, coverage 10X, using the default sequencing substitution, GC% bias, and mutation rate models.

```
python gen_reads.py -r ref.fa -R 101 -o simulated_data
``` 

The most commonly added options are --pe, --bam, --vcf, and -c. 


| Option              | Description                                                                                                                                                                                                                   |
|---------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| @help                 | Displays usage information.
| @reference <str>      | Absolute path to input reference fasta file. | required = yes
| @partition_mode <str> | How to partition the reference for analysis. By default, NEAT will attempt to process one contig per thread. However, if you have very large fasta files, you will see additional runtime benefit from choosing the subdivision method, which will split the contigs up into equal sizes for processing. If you need further speedups and have access to a distributed system you can use a shell script wrapper around NEAT to split the fasta into contigs, then join the results later. NEAT does not feature translocations, so this will not affect NEAT's output. Note that subdivision will only activate for number of threads > 1. | required = no | default = chrom | possible values: chrom, subdivision
| @read_len <int>       | Read length of the reads in the fastq output. Only required if @produce_fastq is set to true. required = no | default = 101
| @threads <int>        | Number of threads to request for NEAT. The recommended amount is the number of chromosomes in your input fasta plus 1. required: no | default = 1
| @output <str>         | Output prefix. Use this option to specify where and what to call output files. | required = yes
| @coverage <float>     | Average coverage for the entire genome. | required: no | default = 10.0
| @error_model <str>    | Absolute path to file with sequencing error model. required: no | default: <NEAT_DIR>/models/errorModel_default.pickle.gz
| @avg_seq_error <float>| Average sequencing error rate for the sequencing machine. Must be between 0.0 and 0.3. | required = no
| @rescale__qualities <bool> | This scales the quality scores to match the desired average sequencing error rate specified by @avg_seq_error. | required = no | default = false
| @ploidy <int>         | Desired ploidy. | required = no | default = 2
| @include_vcf <str>    | Absolute path to vcf file containing variants that will always be included, regardless of genotype and filter. You can pre-filter your vcf for these fields before inputting it if this is not the desired behavior. | required = no
| @target_bed <str>     | Absolute path to bed file containing reference regions that the simulation should target. | required = no
| @discard_bed <str>    | Absolute path to bed file containing reference regions that the simulation should discard. | required = no
| @off_target_coverage <float> | Scalar value for coverage in regions outside the targeted bed. For example, 0.5 would yield roughly half the coverage as the on target areas. | required = no | default = 0.0                                                                     |
| @discard_offtarget <bool> | Whether to discard areas outside the targeted bed region. By default, this is set to false and NEAT will use a different model for off-target regions but still include them in the final output. | required = no | default = false
| @discard_bed <str> | Absolute path to bed file containing reference regions that the simulation should discard. | required = no
| @mutation_model <str> | Absolute path to the mutation model pickle file. Omitting this value will cause NEAT to use the default model, with some standard parameters, and generally uniform biases. | required = no | default: None
| @mutation_rate <float>| Average mutation rate. The mutation rate model is rescaled to make this the average value. Must be between 0 and 0.3.
| @mutation_bed <str>   | Absolute path to a bed file with mutation rates by region. Rates must be in the fourth column and be of the form "mut_rate=x.xx". Rates must be between 0.00 and 0.03. | required = no
| @n_cutoff = <int>	 | Cutoff score for changing a base-call to N.  Below the cutoff, base will be changed to N. If omitted, all bases will be shown, regardless of quality score. | required = no
| @input_vcf <str>         | Input VCF file. Variants from this VCF will be inserted into the simulated sequence with 100% certainty. | required = no
| @paired_ended_data <int> <int> | Paired-end fragment length mean and standard deviation. To produce paired end data, @paired_ended_model must be specified. | required = no
| @paired_ended_model <str>    | Whether the output should be paired ended. For certain conditions (i.e., vcf only or fasta only), this will be ignored. If this is true, then there must be an included fragment length model output from compute_fraglen.py or a mean and standard deviation by declaring values for @fragment_mean and @fragment_std_dev. | required = no | default = false
| @gc_model <str>    | Absolute path to GC content model generated by compute_gc.py. | required = no | default: <NEAT_DIR>/models/gcBias_default.pickle.gz
| @fragment_model <str> | Absolute path to a pickle file containing the fragment length model output from compute_fraglen.py. One such file, fraglenModel_default.pickle.gz is included in the models dir. | required = no
| @fragment_mean <float> | Mean for the paired end fragment length. This only applies if paired-ended is set to true. This number will form the mean for the sample distribution of the fragment lengths in the simulation. Note: This number is REQUIRED if paired_ended is set to true, unless a fragment length model is used. | required = no (unless paired-ended)
| @fragment_st_dev <float> | Standard deviation for the paired end fragment length. This only applies if paired-ended is set to true. This number will form the standard deviation about the mean specified above for the sample distribution of the fragment lengths in the simulation. | required = no (unless paired-ended)
| @produce_bam <bool>  | Whether to produce the golden bam file. This file will contain the reads aligned with the exact region of the genome. | required = no | default = false
| @produce_vcf <bool>  | Whether to produce a vcf file containing all the mutation errors added by NEAT. | required = no | default = false
| @produce_fasta <bool>| Whether to output the mutated fasta. This will output a fasta file with mutations inserted. It does not include sequencing errors or read information. Useful for multigenerational mutations. | required = no | default = false
| @produce_fastq <bool> | Whether to output the fastq(s) of the reads. This is the default output. NEAT will produce 1 fastq for single ended reads or 2 fastqs for paired ended. | required = no | default = true
| @rng <int>          | Set an RNG seed value. Runs using identical RNG values should produce identical results so things like read locations, variant positions, error positions, etc. should be the same. Useful for debugging. | required = no
| @gzip <bool>               | Gzip output FASTQ and VCF. | required = no
| @min_mutations <int> | Set an absolute minimum number of mutations. The program always adds at least 1 mutation. Useful for very small datasets. | required = no
| @discard_offtarget <bool> | Discard reads outside of targeted regions. | required = no
| @rescale_qual      | Rescale quality scores to match @error_model input
| @force_coverage <bool> | If set to true, this will ignore statistical models and force coverage to be constant across the genome. This is considered a debugging feature. | required = no | default = false
| @output_config <bool> | Output config. If yes will output the config file. Useful for repeating runs. | required = no | default = false 
| @debug <bool>    | Turn on debug mode by setting this to true. Debug mode will print certain messages that may help you pinpoint the problem. | required = no | default = false      

## Functionality

![Diagram describing the way that genReads simulates datasets](docs/NEATNEAT.png "Diagram describing the way that genReads simulates datasets")

NEAT produces simulated sequencing datasets. It creates FASTQ files with reads sampled from a provided reference genome, using sequencing error rates and mutation rates learned from real sequencing data. The strength of genReads lies in the ability for the user to customize many sequencing parameters, produce 'golden', true positive datasets, and produce types of data that other simulators cannot (e.g. tumour/normal data).

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
- Parallelized. Can run multiple "partial" simulations in parallel and merge results
- Low memory footprint. Constant (proportional to the size of the reference sequence)

## Examples

The following commands are examples for common types of data to be generated. The simulation uses a reference genome in fasta format to generate reads of 126 bases with default 10X coverage. Outputs paired fastq files, a BAM file and a VCF file. The random variants inserted into the sequence will be present in the VCF and all of the reads will show their proper alignment in the BAM. Unless specified, the simulator will also insert some "sequencing error" -- random variants in some reads that represents false positive results from sequencing.

### Whole genome simulation
Simulate whole genome dataset with random variants inserted according to the default model. 

```
python gen_reads.py                  \
        -r hg19.fa                  \
        -R 126                      \
        -o /home/me/simulated_reads \
        --bam                       \
        --vcf                       \
        --pe 300 30
```

### Targeted region simulation
Simulate a targeted region of a genome, e.g. exome, with -t.

```
python gen_reads.py                  \
        -r hg19.fa                  \
        -R 126                      \
        -o /home/me/simulated_reads \
        --bam                       \
        --vcf                       \
        --pe 300 30                 \
        -t hg19_exome.bed
```

### Insert specific variants
Simulate a whole genome dataset with only the variants in the provided VCF file using -v and -M.

```
python gen_reads.py                  \
        -r hg19.fa                  \
        -R 126                      \
        -o /home/me/simulated_reads \
        --bam                       \
        --vcf                       \
        --pe 300 30                 \
        -v NA12878.vcf              \
        -M 0
```

### Single end reads
Simulate single-end reads by omitting the --pe option.

```
python gen_reads.py                  \
        -r hg19.fa                  \
        -R 126                      \
        -o /home/me/simulated_reads \
        --bam                       \
        --vcf                       
```

### Large single end reads
Simulate PacBio-like reads by providing an error model.

```
python gen_reads.py                             \
	-r hg19.fa                                  \
	-R 5000                                     \
	-e models/errorModel_pacbio_toy.pickle.gz   \
	-E 0.10                                     \
	-o /home/me/simulated_reads        
```

# Utilities	
Several scripts are distributed with gen_reads that are used to generate the models used for simulation.

## compute_gc.py

Computes GC% coverage bias distribution from sample (bedrolls genomecov) data.
Takes .genomecov files produced by BEDtools genomeCov (with -d option).

```
bedtools genomecov
        -d                          \
        -ibam normal.bam            \
        -g reference.fa
```

```
python compute_gc.py                \
        -r reference.fa             \
        -i genomecovfile            \
        -w [sliding window length]  \
        -o /path/to/prefix
```

## compute_fraglen.py

Computes empirical fragment length distribution from sample data.
Takes SAM file via stdin:

    python computeFraglen.py \
        -i input.bam         \
        -o /prefix/for/output

and creates fraglen.pickle.gz model in working directory.

## generate_mutation_model.py

Takes references genome and VCF file to generate mutation models:

```
python gen_mut_model.py               \
        -r hg19.fa                  \
        -m inputVariants.vcf        \
        -o /home/me/models
```

Trinucleotides are identified in the reference genome and the variant file. Frequencies of each trinucleotide transition are calculated and output as a pickle (.p) file.

| Option          | Description                                                                  |
|-----------------|------------------------------------------------------------------------------|
| -r <str>        | Reference file for organism in FASTA format. Required                        |
| -m <str>        | Mutation file for organism in VCF format. Required                           |
| -o <str>        | Path to output file and prefix. Required.                                    |
| --bed           | Flag that indicates you are using a bed-restricted vcf and fasta (see below) |
| --save-trinuc   | Save trinucleotide counts for reference                                      |
| --human-sample  | Use to skip unnumbered scaffolds in human references                         |
| --skip-common   | Do not save common snps or high mutation areas                               |

Note that if you have a bed input, you will need to have bedtools installed in your environment, in addition to the 
pybedtools python package. We recommend using Anaconda:

```
conda install -c bioconda bedtools
```

See the [bedtools documentation](https://bedtools.readthedocs.io/en/latest/) for more information.

## genSeqErrorModel.py

Generates sequence error model for gen_reads.py -e option.
This script needs revision, to improve the quality-score model eventually, and to include code to learn sequencing errors from pileup data.

```
python genSeqErrorModel.py                            \
        -i input_read1.fq (.gz) / input_read1.sam     \
        -o /output/prefix                             \
        -i2 input_read2.fq (.gz) / input_read2.sam    \
        -p input_alignment.pileup                     \
        -q quality score offset [33]                  \
        -Q maximum quality score [41]                 \
        -n maximum number of reads to process [all]   \
        -s number of simulation iterations [1000000]  \
        --plot perform some optional plotting
```

## plotMutModel.py

Performs plotting and comparison of mutation models generated from genMutModel.py.

```
python plotMutModel.py                                                  \
        -i model1.pickle.gz [model2.pickle.gz] [model3.pickle.gz]...    \
        -l legend_label1 [legend_label2] [legend_label3]...             \
        -o path/to/pdf_plot_prefix
```

## vcf_compare_OLD.py

Tool for comparing VCF files.

```
python vcf_compare_OLD.py
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
ICGC's "Access Controlled Data" documention can be found at http://docs.icgc.org/access-controlled-data. To have access to controlled germline data, a DACO must be
submitted. Open tier data can be obtained without a DACO, but germline alleles that do not match the reference genome are masked and replaced with the reference
allele. Controlled data includes unmasked germline alleles.




