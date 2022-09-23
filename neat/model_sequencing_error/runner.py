"""
Creates a model of sequencing errors in the dataset. One assumption of this model is that low quality scores indicate
sequencing errors. We should verify this is true in practice and/or how that assumption lines up with variant callers.

I wonder if it would be possible to get a list of all items a variant caller deems errors.
"""

import gzip
import pickle
import numpy as np
import logging

from pathlib import Path

from .utils import take_closest, parse_file
from ..common import validate_output_path

__all__ = [
    "model_seq_err_runner"
]

_LOG = logging.getLogger(__name__)


def model_seq_err_runner(
        file: str | Path,
        offset: int,
        qual_scores: list | int,
        max_reads: int,
        num_iters: int,
        plot: bool,
        overwrite: bool,
        output_prefix: str
):
    """
    Sets up and calls the sequencing error modeling core code.

    :param file: This is the input file, which can be in fastq (gzipped or not), sam, or bam format.
    :param offset: This is the quality offset. Illumina machines use the sanger model qual score + 33 to
        derive the character for the quality array.
    :param qual_scores: These are the possible quality scores for the dataset. This can either be a list or a single
        integer.
    :param max_reads: The maximum number or reads to process. This can speed up the time taken to create the model,
        at the expense of accuracy.
    :param num_iters: The number of simulation iterations to run. This is to model the probability.
    :param plot: run optional plotting.
    :param overwrite: True to overwrite input, mainly a debugging option
    :param output_prefix: The name of the file to write the output.
    """

    input_file = Path(file)
    gzipped = False
    common_fastq_suffixes = ['.fq', '.fastq']
    if ".bam" in input_file.suffixes:
        filetype = 'bam'
    elif any([x for x in input_file.suffixes if x in common_fastq_suffixes]):
        filetype = "fastq"
        if ".gz" in input_file.suffixes:
            gzipped = True
    elif '.sam' in input_file.suffixes:
        filetype = 'sam'
    else:
        raise ValueError("Unknown filetype. Please enter a fastq bam/sam file, with common "
                         "suffixes [.sam, .bam, .fq(.gz) or .fastq(.gz)] (gzipping is optional).")

    _LOG.info(f'Input File: {str(input_file)}')
    _LOG.info(f'Filetype: {filetype}')
    _LOG.info(f'File is gzipped: {gzipped}')

    message = f"[NOT 33? You animal!] {offset}" if offset != 33 else "33"
    _LOG.info(f"Quality offset: {message}")

    final_quality_scores: list
    if len(qual_scores) == 1:
        final_quality_scores = list(range(qual_scores[0] + 1))
    else:
        final_quality_scores = qual_scores

    _LOG.info(f'Quality scores: {final_quality_scores}')
    num_records_to_process = "all" if max_reads == -1 else str(max_reads)
    _LOG.info(f'Maximum number of records to process: {num_records_to_process}')
    _LOG.info(f'Number of simulation iterations to run: {num_iters}')
    _LOG.info(f"Plot the data? {plot}")
    _LOG.info(f'Overwrite existing data? {overwrite}')
    output_file = Path(output_prefix).with_suffix(".p.gz")
    _LOG.info(f'Writing output to: {output_file}')

    parse_file(input_file, final_quality_scores, offset, max_reads, num_iters)

    _LOG.info("All done, have a nice day!")


    # (infile, outfile, off_q, q_scores, max_reads, n_samp) = (args.i, args.o, args.q, args.Q, args.n, args.s)
    # # (infile2, pile_up) = (args.i2, args.p)
    # infile2 = None
    # pile_up = None
    #
    # plot_stuff = args.plot
    #
    # # q_scores = range(real_q)
    # if infile2 is None:
    #     (init_q, prob_q, avg_err) = parse_file(infile, q_scores, off_q, max_reads, n_samp)
    # else:
    #     (init_q, prob_q, avg_err1) = parse_file(infile, q_scores, off_q, max_reads, n_samp)
    #     (init_q2, prob_q2, avg_err2) = parse_file(infile2, q_scores, off_q, max_reads, n_samp)
    #     avg_err = (avg_err1 + avg_err2) / 2.
    #
    # # Embed some default sequencing error parameters
    # # sequencing substitution transition probabilities
    # sse_prob = [[0., 0.4918, 0.3377, 0.1705],
    #             [0.5238, 0., 0.2661, 0.2101],
    #             [0.3754, 0.2355, 0., 0.3890],
    #             [0.2505, 0.2552, 0.4942, 0.]]
    # # if a sequencing error occurs, what are the odds it's an indel?
    # sie_rate = 0.01
    # # sequencing indel error length distribution
    # sie_prob = [0.999, 0.001]
    # sie_val = [1, 2]
    # # if a sequencing indel error occurs, what are the odds it's an insertion as opposed to a deletion?
    # sie_ins_freq = 0.4
    # # if a sequencing insertion error occurs, what's the probability of it being an A, C, G, T...
    # sie_ins_nucl = [0.25, 0.25, 0.25, 0.25]
    #
    # # Otherwise we need to parse a pileup and compute statistics!
    # if pile_up:
    #     print('\nPileup parsing coming soon!\n')
    #     exit(1)
    # else:
    #     print('Using default sequencing error parameters...')
    #
    # err_params = [sse_prob, sie_rate, sie_prob, sie_val, sie_ins_freq, sie_ins_nucl]
    #
    # # finally, let's save our output model
    # outfile = pathlib.Path(outfile).with_suffix(".p")
    # print('saving model...')
    # if infile2 is None:
    #     # pickle.dump({quality_scores:[], quality_scores_probabilities:[]})
    #     pickle.dump([init_q, prob_q, q_scores, off_q, avg_err, err_params], open(outfile, 'wb'))
    # else:
    #     pickle.dump([init_q, prob_q, init_q2, prob_q2, q_scores, off_q, avg_err, err_params], open(outfile, 'wb'))

