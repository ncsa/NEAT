#!/usr/bin/env python

"""
The purpose of this script is to run gen reads in multithreaded mode. The way it will work is to
split the reference fasta up into chromosomes and then run gen_reads in one thread per chromosome mode.
Finally, it will recombine the final output into one.

Note that I'm planninng to hijack NEAT's features for this.
"""

import argparse
import bisect
import gzip
import logging
import pathlib
import random
import sys
import time
import tempfile
import pandas as pd
import numpy as np
import glob
from heapq import merge
import os
import atexit

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from source.Models import Models
from source.Options import Options
from source.bed_func import parse_bed
from source.constants_and_defaults import ALLOWED_NUCL
from source.constants_and_defaults import VERSION
from source.error_handling import premature_exit, log_mssg
from source.output_file_writer import OutputFileWriter
from source.run_neat import generate_variants, generate_fasta, generate_reads, generate_bam
from source.vcf_func import parse_input_vcf


class StdoutFilter(logging.Filter):
    # Setup filter for logging options
    def filter(self, record):
        return record.levelno == logging.INFO


# useful functions for the rest
def print_start_info():
    """
    Prints a NEAT start message with version number and logs the start time

    :return: logs the time this function was called as start time.
    """
    starttime = time.time()
    print('\n-----------------------------------------------------------')
    log_mssg(f"NEAT multithreaded version, v{VERSION}, is running.", 'info')
    log_mssg(f'Started: {str(starttime)}.', 'info')
    return starttime


def print_end_info(output_class, cancer_output_class, starting_time):
    endtime = time.time()
    list_to_iterate = output_class.files_to_write
    if cancer_output_class:
        list_to_iterate += cancer_output_class.files_to_write
    log_mssg("Wrote the following files:", 'info')
    i = 1
    for item in list_to_iterate:
        log_mssg(f"\t{i}. {item}", 'info')
        i += 1
    log_mssg("NEAT finished successfully.", "info")
    log_mssg(f"Total runtime: {str(endtime - starting_time)}", 'info')
    print("-------------------------------------------------------------------------")


def print_configuration(args, options):
    """
    Combines the relevant parts of the input args and the options file to print out a
    list of the configuration parameters. Useful for reproducibility.
    """
    log_mssg(f'Run Configuration...', 'INFO')
    potential_filetypes = ['vcf', 'bam', 'fasta', 'fastq']
    log = ''
    for suffix in potential_filetypes:
        key = f'produce_{suffix}'
        if options.args[key]:
            log += f'{suffix} '
    log_mssg(f'Producing the following files: {log.strip()}', 'INFO')
    log_mssg(f'Input file: {options.reference}', 'INFO')
    log_mssg(f'Output files: {args.output}.<{log}>', 'INFO')
    if options.threads == 1:
        log_mssg(f"Single threading - 1 thread.", 'info')
    else:
        log_mssg(f'Multithreading coming soon!!', 'info')
        log_mssg(f"Single threading - 1 thread.", 'info')
        # We'll work on mulitthreading later...
        # print_and_log(f'Multithreading - {options.threads} threads', 'info')
    if options.paired_ended:
        log_mssg(f'Running in paired-ended mode.', 'INFO')
        if options.fragment_model:
            log_mssg(f'Using fragment length model: {options.fragment_model}', 'INFO')
        else:
            log_mssg(f'Using fragment model based on mean={options.fragment_mean}, '
                     f'st dev={options.fragment_st_dev}', 'INFO')
    else:
        log_mssg(f'Running in single-ended mode.', 'INFO')
    log_mssg(f'Using a read length of {options.read_len}', 'INFO')
    log_mssg(f'Average coverage: {options.coverage}', 'INFO')
    log_mssg(f'Using error model: {options.error_model}', 'INFO')
    if options.avg_seq_error:
        log_mssg(f'User defined average sequencing error rate: {options.avg_seq_error}.', 'INFO')
    if options.rescale_qualities:
        log_mssg(f'Quality scores will be rescaled to match avg seq error rate.', 'INFO')
    log_mssg(f'Ploidy value: {options.ploidy}', 'INFO')
    if options.include_vcf:
        log_mssg(f'Vcf of variants to include: {options.include_vcf}', 'INFO')
    if options.target_bed:
        log_mssg(f'BED of regions to target: {options.target_bed}', 'INFO')
        log_mssg(f'Off-target coverage rate: {options.off_target_coverage}', 'INFO')
        log_mssg(f'Discarding off target regions: {options.discard_offtarget}', 'INFO')
    if options.discard_bed:
        log_mssg(f'BED of regions to discard: {options.discard_bed}', 'INFO')
    if options.mutation_model:
        log_mssg(f'Using mutation model in file: {options.mutation_model}', 'INFO')
    if options.mutation_rate:
        log_mssg(f'Rescaling average mutation rate to: {options.mutation_rate}', 'INFO')
    if options.mutation_bed:
        log_mssg(f'BED of mutation rates of different regions: {options.mutation_bed}', 'INFO')
    if options.n_cutoff:
        log_mssg(f'N-cutoff quality score: {options.n_cutoff}', 'INFO')
    if options.gc_model:
        log_mssg(f'Using GC model: {options.gc_model}', 'INFO')
    if options.force_coverage:
        log_mssg(f'Ignoring models and forcing coverage value.', 'INFO')
    log_mssg(f'Debug Mode Activated.', 'debug')
    log_mssg(f'RNG seed value: {options.rng_value}', 'debug')


def find_file_breaks(threads: int, mode: str, reference_index: dict) -> dict:
    """
    Returns a dictionary with the chromosomes as keys
    For the chrom method, the value for each key will just  be "all"
    whereas for subdivison, the value for each key should be a list of indices that partition the
    sequence into roughly equal sizes.
    :param threads: number of threads for this run
    :param mode: partition mode for this run (chrom or subdivision)
    :param reference_index: a dictionary with chromosome keys and sequence values
    :return: a dictionary containing the chromosomes as keys and either "all" for valuse, or a list of indices

    >>> index = {'chr1': "ACCATACACGGGCAACACACGTACACATTATACC"}
    >>> find_file_breaks(5, "subdivision", index)
    {'chr1': [(0, 6), (6, 12), (12, 18), (18, 24), (24, 34)]}
    >>> find_file_breaks(4, 'subdivision', index)
    {'chr1': [(0, 8), (8, 16), (16, 24), (24, 34)]}
    >>> find_file_breaks(5, 'chrom', index)
    {'chr1': [(0, 34)]}
    >>> find_file_breaks(2, 'subdivision', index)
    {'chr1': [(0, 17), (17, 34)]}
    """
    partitions = {}
    if mode.lower() == "chrom" or threads == 1:
        for contig in reference_index.keys():
            partitions[contig] = [(0, len(reference_index[contig]))]
    elif mode.lower() == "subdivision":
        # Add items one at a time to partition list until the total length is greater than delta.
        for contig in reference_index:
            contig_length = len(reference_index[contig])
            delta = contig_length // threads

            if contig not in partitions:
                partitions[contig] = []

            breakpoints = list(range(0, contig_length, delta))

            # Since we know the first index will be zero, we can skip the first item in the breakpoints list
            # And since we want the last partition to grab the rest, we'll stop short and add it manually
            for index in breakpoints:
                if index + delta <= contig_length:
                    partitions[contig].append((index, index + delta))
                else:
                    # Have to extend the last one so we don't get a tiny read we can't process
                    partitions[contig][-1] = (partitions[contig][-1][0], contig_length)

        log_mssg(f'breaks = {partitions}', 'debug')
    else:
        log_mssg("Invalid partition mode. Must be either chrom or subdivision.", 'error')
        premature_exit(1)
    return partitions


def main(raw_args=None):
    """
    Parallel main function. Takes args and parses the ones needed to start multiprocessing. The rest will get passed
    along to the main processing part of gen_reads.
    """

    parser = argparse.ArgumentParser(description=f'This script runs gen_reads v{VERSION}',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, )
    parser.add_argument('-c', '--config', default='neat_config.txt', dest='conf', metavar='config',
                        help="Any extra arguments valid for NEAT")
    parser.add_argument('-o', '--output', required=True,
                        help="Prefix for the output. Can include a path before it.")

    args = parser.parse_args(raw_args)

    # This will by default look for a file called "neat_config.txt" in the current working dir
    args.conf = pathlib.Path(args.conf)

    log_file = f'{args.output}.log'
    if pathlib.Path(log_file).exists():
        os.remove(log_file)

    starttime = print_start_info()

    if not args.conf.is_file():
        log_mssg(f'Configuration file ({args.conf}) cannot be found.', 'error')
        premature_exit(1)

    # Reads in the user-entered options from the config file and performs
    # some basic checks and corrections.
    options = Options(args.conf)

    if not options.debug:
        file_handler.setLevel(logging.INFO)

    print_configuration(args, options)

    """
    Model preparation
    
    Read input models or default models, as specified by user.
    """
    log_mssg("Reading Models...", 'info')

    models = Models(options)

    """
    Process Inputs
    """
    log_mssg("Processing inputs...", 'info')

    log_mssg(f'Reading {options.reference}...', 'info')

    reference_index = SeqIO.index(str(options.reference), 'fasta')
    log_mssg(f'Reference file indexed.', 'debug')

    # there is not a reference_chromosome in the options defs because
    # there is nothing that really needs checking, so this command will just
    # set the value. This might be an input in the future, for example if you want
    # to only simulate certain chroms.
    options.set_value('reference_chromosomes', list(reference_index.keys()))

    # TODO maybe take this out. I feel like we're overthinking chr, but I could be wrong
    # this parameter will ensure we account for 'chr' in the output.
    # It is True if they all start with "chr" and False otherwise.
    begins_with_chr = all(k.startswith('chr') for k in options.reference_chromosomes)

    sample_names = []
    input_variants = {}
    if options.include_vcf:
        log_mssg(f"Reading input VCF...", 'info')
        if options.cancer:
            (sample_names, input_variants) = parse_input_vcf(options.include_vcf,
                                                             options.ploidy,
                                                             models.mutation_model['homozygous_freq'],
                                                             reference_index,
                                                             tumor_normal=True)
            # TODO figure out what these were going to be used for
            tumor_ind = sample_names.index('tumor_sample')
            normal_ind = sample_names.index('normal_sample')
        else:
            (sample_names, input_variants) = parse_input_vcf(options.include_vcf,
                                                             options.ploidy,
                                                             models.mutation_model['homozygous_freq'],
                                                             reference_index)

        log_mssg("Finished reading @include_vcf file.", "info")

    # parse input targeted regions, if present.
    if options.target_bed or options.discard_bed or options.mutation_bed:
        log_mssg(f"Reading input bed files.", 'info')

    # Note if any bed is empty, parse_bed just returns a dict of chromosomes with empty list values

    target_regions_dict = parse_bed(options.target_bed, options.reference_chromosomes,
                                    begins_with_chr, False)

    discard_regions_dict = parse_bed(options.discard_bed, options.reference_chromosomes,
                                     begins_with_chr, False)

    mutation_rate_dict = parse_bed(options.mutation_bed, options.reference_chromosomes,
                                   begins_with_chr, True)

    log_mssg(f'Finished reading input beds.', 'debug')

    """
    Initialize Output Files
    """
    log_mssg("Initializing output files...", 'info')

    # Prepare headers
    bam_header = None
    vcf_header = None
    if options.produce_bam:
        bam_header = reference_index
    if options.produce_vcf:
        vcf_header = [options.reference]

    out_prefix_name = pathlib.Path(args.output).resolve().name
    out_prefix_parent_dir = pathlib.Path(pathlib.Path(args.output).resolve().parent)
    if not out_prefix_parent_dir.is_dir():
        log_mssg(f'Creating output dir: {out_prefix_parent_dir}', 'debug')
        out_prefix_parent_dir.mkdir(parents=True, exist_ok=True)

    # Creates files and sets up objects for files that can be written to as needed.
    # Also creates headers for bam and vcf.
    # We'll also keep track here of what files we are producing.
    if options.cancer:
        output_normal = out_prefix_parent_dir / f'{out_prefix_name}_normal'
        output_tumor = out_prefix_parent_dir / f'{out_prefix_name}_tumor'
        output_file_writer = OutputFileWriter(output_normal,
                                              bam_header=bam_header,
                                              vcf_header=vcf_header,
                                              options_file=options
                                              )
        output_file_writer_cancer = OutputFileWriter(output_tumor,
                                                     bam_header=bam_header,
                                                     vcf_header=vcf_header,
                                                     options_file=options
                                                     )
    else:
        outfile = out_prefix_parent_dir / out_prefix_name
        output_file_writer = OutputFileWriter(outfile,
                                              bam_header=bam_header,
                                              vcf_header=vcf_header,
                                              options_file=options
                                              )
        output_file_writer_cancer = None

    log_mssg(f'Output files ready for writing.', 'debug')

    """
    Begin Analysis
    """
    log_mssg("Beginning analysis...", 'info')

    # Find break points in the input file.
    # TODO Debug:
    breaks = find_file_breaks(options.threads, options.partition_mode, reference_index)

    if not breaks:
        # Printing out summary information and end time
        log_mssg("Found no chromosomes in reference.", 'debug')
        print_end_info(output_file_writer, output_file_writer_cancer, starttime)
        sys.exit(0)

    log_mssg("Input file partitioned.", 'debug')

    """
    REDO with mpire. The package can map a function to an iterable, so we will use the breaks as the iterable
    And recast neat as a function (run_neat) that will take in the items from that iterable and process them.
    I think the core of run_neat should be taking in a sequence and returning a mutated sequence. The core of NEAT
    at the moment tries to be everything to everyone.
    """
    log_mssg("Beginning simulation", 'info')
    # these will be the features common to each contig, for multiprocessing
    common_features = {}

    out_prefix = f'{out_prefix_parent_dir}/{out_prefix_name}'

    vcf_files = []
    fasta_records = []
    fastq_files = []
    bam_files = []

    temporary_dir = tempfile.TemporaryDirectory()

    for contig in breaks:
        log_mssg(f'Mutating: {contig}...', 'info')

        contig_variants = {x[1]: input_variants[x] for x in input_variants if x[0] == contig}
        reference = reference_index[contig]

        # Since we're only running single threaded for now:
        threadidx = 1

        tmp_dir_path = pathlib.Path(temporary_dir.name).resolve()

        # We need the temp vcf file no matter what
        temp_vcf_filename = tmp_dir_path / f"{out_prefix_name}_tmp_{contig}_{threadidx}.vcf"
        log_mssg(f'temp_vcf_filename = {temp_vcf_filename}', 'debug')

        if threadidx == 1:
            # init_progress_info()
            pass

        chrom_vcf_file, tmp_vcf_file = generate_variants(reference,
                                                         contig,
                                                         temp_vcf_filename,
                                                         target_regions_dict[contig],
                                                         discard_regions_dict[contig],
                                                         mutation_rate_dict[contig],
                                                         contig_variants,
                                                         models,
                                                         options,
                                                         out_prefix)

        vcf_files.append(chrom_vcf_file)

        if options.produce_fastq:
            # Note that generate_reads returns None for r2 if single ended
            chrom_fastq_r1_file, chrom_fastq_r2_file = generate_reads(reference,
                                                                      models,
                                                                      tmp_vcf_file,
                                                                      tmp_dir_path,
                                                                      options,
                                                                      contig)
            fastq_files.append([chrom_fastq_r1_file, chrom_fastq_r2_file])

        if options.produce_fasta:
            chrom_fasta_dict = generate_fasta(reference, options, tmp_vcf_file,
                                              tmp_dir_path, contig)
            fasta_records.append(SeqRecord(Seq(chrom_fasta_dict), id=reference.id, description=reference.description))

        if options.produce_bam:
            chrom_bam_file = generate_bam(reference, options, tmp_vcf_file, tmp_dir_path, contig)
            bam_files.append(chrom_bam_file)

    if options.produce_vcf:
        log_mssg("Outputting complete golden vcf.", 'info')
        chunks = []
        # TODO double check that breaks is in the same order as the input fasta
        for vcf in vcf_files:
            chunks += [gzip.open(vcf, 'r')]

        with gzip.open(output_file_writer.vcf_fn, 'ab') as vcf_out:
            vcf_out.writelines(merge(*chunks))

        for item in chunks:
            item.close()
            os.remove(item.name)

    if options.produce_fasta:
        log_mssg("Outputting final fasta file", 'info')

        with gzip.open(output_file_writer.fasta_fn, 'at') as fasta_out:
            SeqIO.write(fasta_records, fasta_out, 'fasta')


        # for fasta in fasta_files:
        #     in_file = open(fasta, 'r')
        #     fasta_out.write(in_file.read())
        #     in_file.close()
        #
        # fasta_out.close()

    # End info
    print_end_info(output_file_writer, output_file_writer_cancer, starttime)


if __name__ == '__main__':

    now = time.localtime()
    now = time.strftime('%Y_%m_%d_%H%M', now)

    log_dir = os.getcwd()

    log_name = f'{log_dir}/{now}_NEAT.log'

    neat_log = logging.getLogger()
    neat_log.setLevel(logging.DEBUG)
    neat_log.handlers = []

    formatter = logging.Formatter('%(asctime)s %(levelname)s %(threadName)s: %(message)s', datefmt='%Y/%m/%d %H:%M')

    stdout_handler = logging.StreamHandler(stream=sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    stdout_handler.addFilter(StdoutFilter())
    stdout_handler.setFormatter(formatter)

    neat_log.addHandler(stdout_handler)

    stderr_handler = logging.StreamHandler(stream=sys.stderr)
    stderr_handler.setFormatter(formatter)
    stderr_handler.setLevel(logging.WARNING)

    neat_log.addHandler(stderr_handler)

    file_handler = logging.FileHandler(log_name)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    neat_log.addHandler(file_handler)

    main()
