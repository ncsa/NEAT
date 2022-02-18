import time
from mpire import WorkerPool
from source.Options import Options
import numpy as np
import pandas as pd
import networkx as nx


def close_temp_files(tmp_files):
    for file_handle in tmp_files:
        file_handle.close()


def close_temp_dir(temp_dir):
    temp_dir.cleanup()


class Job():
    def __init__(self, ref_idx, out_prefix_name, target_regions, discard_regions,
                 mutation_rate_regions, input_variants, models, options):
        self.ref_idx = ref_idx
        self.out_prefix_name = out_prefix_name
        self.target_regions = target_regions
        self.discard_regions = discard_regions
        self.mutation_rate_regions = mutation_rate_regions
        self.input_variants = input_variants
        self.models = models
        self.options = options


def run_neat(params, job_start, job_end):
    """
    >>> job = (0,10)
    >>> common = Job({}, 'hello', None, None, None, None, None, None)
    >>> x = run_neat(common, job[0], job[1])
    hello
    >>> print(x)
    Hello World
    """
    print(params.out_prefix_name)
    for index in range(job_start, job_end):
        time.sleep(0.1)
    return "Hello World"


def execute_neat(options, jobs, params):
    """
    >>> job = ((10, 15), (15, 20), (20, 25))
    >>> common = Job({'chr1': "ACCTAGACA"}, 'hello', pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), None, None)
    >>> options = Options()
    >>> options.threads = 3
    >>> execute_neat(options, job, common)
    ['Hello World', 'Hello World', 'Hello World']
    All done, have a nice day!
    """
    with WorkerPool(n_jobs=options.threads, shared_objects=params) as pool:
        results = pool.map(run_neat, jobs, iterable_len=len(jobs), progress_bar=True)

    print(results)
    print("All done, have a nice day!")
    # TODO idea what if we loop over the chromosome once and add mutations. Then do sequence errors during the reads
    # TODO add structural variants here.

    # Step 1: Create a VCF of variants (mutation and sequencing error variants)
    # Step 2 (optional): create a fasta with those variants inserted
    # Step 3 (optional): create a fastq from the mutated fasta
    # Step 4 (optional): Create a golden bam with the reads aligned to the original reference

    # reference = ref_index[chromosome]
    #
    # # These tasks filter the inputs down to relevant areas
    # partitioned_reference = {}
    # if len(partition) == 1:
    #     # In this case, we just use the whole sequence
    #     partitioned_reference[(0, len(ref_index[chromosome]))] = ref_index[chromosome].seq
    # else:
    #     for i in range(len(partition) - 1):
    #         # fetch the coordinates
    #         start = partition[i]
    #         stop = partition[i + 1]
    #         # index is a tuple: (chromosome (string), region start (int), region end (int))
    #         partitioned_reference[(start, stop)] = ref_index[chromosome][start:stop]
    #     # Grab the rest of the sequence as the final partition
    #     partitioned_reference[(partition[-1], len(ref_index[chromosome]))] = \
    #         ref_index[chromosome][partition[-1]:]
    #
    # n_regions = {}
    #
    # input_variants = pd.DataFrame()
    # if not input_variants.empty:
    #     input_variants = input_variants[input_variants.index.isin([chromosome])]
    #
    # target_regions = pd.DataFrame()
    # if not target_regions.empty:
    #     target_regions = target_regions[target_regions.index.isin([chromosome])]
    #
    # discard_regions = pd.DataFrame()
    # if not discard_regions.empty:
    #     discard_regions = discard_regions[discard_regions.index.isin([chrom])]
    #
    # mutation_rate_regions = pd.DataFrame()
    # if not mutation_rate_regions.empty:
    #     mutation_rate_regions = mutation_rate_regions[mutation_rate_regions.index.isin([chromosome])]
    #
    # # TODO Seq error for each group
    # sequencing_error_class = SequencingErrors(options.read_len, models.sequencing_error_model,
    #                                           options.avg_seq_error, options.rescale_qualities,
    #                                           options.debug)
    #
    # # TODO check into multiprocessing.Pool()
    #
    # # Setting up temp files to write to
    # tmp_fasta_fn = None
    # tmp_fastq1_fn = None
    # tmp_fastq2_fn = None
    # tmp_sam_fn = None
    # tmp_vcf_fn = None
    #
    # temporary_dir = tempfile.TemporaryDirectory()
    # tmp_dir_path = pathlib.Path(temporary_dir.name).resolve()
    # if options.produce_bam:
    #     tmp_sam_fn = tmp_dir_path / f"{out_prefix_name}_tmp_records_{threadidx}.tsam"
    #     if options.debug:
    #         print_and_log(f'tmp_sam_fn = {tmp_sam_fn}', 'debug')
    # if options.produce_vcf:
    #     tmp_vcf_fn = tmp_dir_path / f"{out_prefix_name}_tmp_{threadidx}.vcf"
    #     if options.debug:
    #         print_and_log(f'tmp_vcf_fn = {tmp_vcf_fn}', 'debug')
    # if options.produce_fasta:
    #     tmp_fasta_fn = tmp_dir_path / f"{out_prefix_name}_tmp_{threadidx}.fasta"
    #     if options.debug:
    #         print_and_log(f'tmp_fasta_fn = {tmp_fasta_fn}', 'debug')
    # if options.produce_fastq:
    #     tmp_fastq1_fn = tmp_dir_path / f"{out_prefix_name}_tmp_{threadidx}_read1.fq"
    #     if options.debug:
    #         print_and_log(f'tmp_fastq1_fn = {tmp_fastq1_fn}', 'debug')
    #     if options.paired_ended:
    #         tmp_fastq2_fn = tmp_dir_path / f"{out_prefix_name}_tmp_{threadidx}_read2.fq"
    #         if options.debug:
    #             print_and_log(f'tmp_fastq2_fn  = {tmp_fastq2_fn}', 'debug')
    #
    # # Setting up temp files for writing
    # tmp_fasta_outfile = None
    # tmp_fastq1_outfile = None
    # tmp_fastq2_outfile = None
    # tmp_sam_outfile = None
    # tmp_vcf_outfile = None
    # all_tmp_files = []
    #
    # if tmp_fasta_fn:
    #     tmp_fasta_outfile = open(tmp_fasta_fn, 'w')
    #     all_tmp_files.append(tmp_fasta_outfile)
    # if tmp_fastq1_fn:
    #     tmp_fastq1_outfile = open(tmp_fastq1_fn, 'w')
    #     all_tmp_files.append(tmp_fastq1_outfile)
    # if tmp_fastq2_fn:
    #     tmp_fastq2_outfile = open(tmp_fastq2_fn, 'w')
    #     all_tmp_files.append(tmp_fastq2_outfile)
    # if tmp_sam_fn:
    #     tmp_sam_outfile = open(tmp_sam_fn, 'w')
    #     all_tmp_files.append(tmp_sam_outfile)
    # if tmp_vcf_fn:
    #     tmp_vcf_outfile = open(tmp_vcf_fn, 'w')
    #     all_tmp_files.append(tmp_vcf_outfile)
    #
    # if threadidx == 1:
    #     # init_progress_info()
    #     pass
    #
    # def quick_mutate(dna_string: str, length: int):
    #     mutated_sequence = ""
    #     quality_string = ""
    #     range_to_iterate = min(length, len(dna_string))
    #     for i in range(range_to_iterate):
    #         if random.random() < 0.01:
    #             mutated_sequence += random.choice(ALLOWED_NUCL)
    #         else:
    #             mutated_sequence += dna_string[i]
    #         quality_string += chr(random.randint(2, 40) + 33)
    #     return mutated_sequence, quality_string
    #
    # for i in range(len(partition)):
    #     # This is filler code until we do the actual processing.
    #     tlen = 300
    #     pos = 1
    #     for read_num in range(1000):
    #         qname = f'{out_prefix_name}-{chromosome}-{read_num}'
    #         flag = 0
    #         rname = chromosome
    #         mapq = 70
    #         rnext = '='
    #         pnext = 1
    #         # SAM is 1-based for annoying reasons, therefore we have to subtract 1 to get the correct positions
    #         reference_sequence = \
    #             partitioned_reference[list(partitioned_reference.keys())[0]].seq[pos - 1:pos + tlen]
    #         seq, qual = quick_mutate(reference_sequence, 101)
    #         line_to_write = f'{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{reference_sequence}' \
    #                         f'\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\n'
    #         tmp_sam_outfile.write(line_to_write)
    #         pos += 1
    # tmp_sam_outfile.close()
    # shutil.copy(tmp_sam_fn, f'/home/joshfactorial/Documents/temp_{threadidx}.tsam')
    # tmp_sam_outfile.close()