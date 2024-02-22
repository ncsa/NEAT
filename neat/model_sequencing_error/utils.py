"""
Utilities to generate the sequencing error model
"""

import logging
import numpy as np
# TODO implement plotting
# import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
from scipy.stats import mode
from ..common import open_input
from ..models import take_closest

__all__ = [
    "parse_file"
]

_LOG = logging.getLogger(__name__)


def convert_quality_string(qual_str: str, offset: int):
    """
    Converts a plain quality string to a list of numerical equivalents

    :param qual_str: The string to convert
    :param offset: the quality offset for conversion for this fastq
    :return list: a list of numeric quality scores
    """
    ret_list = []
    for i in range(len(qual_str)):
        try:
            ret_list.append(ord(qual_str[i]) - offset)
        except ValueError:
            raise ValueError("improperly formatted fastq file")

    return ret_list


def expand_counts(count_array: list, scores: list):
    """
    Expands a counting list out into the full tally

    :param count_array: the list to expand
    :param scores: The factors by which to expand the list
    :return np.ndarray: a one-dimensional array reflecting the expanded count
    """
    if len(count_array) != len(scores):
        raise ValueError("Count array and scores have different lengths.")

    ret_list = []
    for i in range(len(count_array)):
        ret_list.extend([scores[i]] * count_array[i])

    return np.array(ret_list)


def parse_file(input_file: str, quality_scores: list, max_reads: int, qual_offset: int, readlen: int):
    """
    Parses an individual file for statistics

    :param input_file: The input file to process
    :param quality_scores: A list of potential quality scores
    :param max_reads: Max number of reads to process for this file
    :param qual_offset: The offset score for this fastq file. We assume the Illumina default of 33.
    :param readlen: The read length for these datasets. If 0, then we will determine it.
    :return:
    """

    _LOG.info(f'reading {input_file}')

    if not readlen:
        readlens = []

        # takes too long and uses too much memory to read all of them, so let's just get a sample.
        with open_input(input_file) as fq_in:
            i = 0
            # Count the first 1000 lines
            while i < 1000:
                i += 1
                for _ in (0, 1, 2):
                    fq_in.readline()
                line = fq_in.readline().strip()
                readlens.append(len(line))

        readlens = np.array(readlens)

        # Using the statistical mode seems like the right approach here. We expect the readlens to be roughly the same.
        readlen_mode = mode(readlens, axis=None, keepdims=False)
        if readlen_mode.count < (0.5 * len(readlens)):
            _LOG.warning("Highly variable read lengths detected. Results may be less than ideal.")
        if readlen_mode.count < 20:
            raise ValueError(f"Dataset is too scarce or inconsistent to make a model. Try a different input.")
        read_length = int(readlen_mode.mode)

    else:
        read_length = readlen

    _LOG.debug(f'Read len of {read_length}, over {1000} samples')

    _LOG.info(f"Reading {max_reads} records...")
    temp_q_count = np.zeros((read_length, len(quality_scores)), dtype=int)
    qual_score_counter = {x: 0 for x in quality_scores}
    # shape_curves = []
    quarters = max_reads//4

    records_read = 0
    wrong_len = 0
    end_of_file = False
    # SeqIO eats up way too much memory for larger fastqs, so we're trying to read the file in line by line here
    with open_input(input_file) as fq_in:
        while records_read < max_reads:

            # We throw away 3 lines and read the 4th, because that's fastq format
            for _ in (0, 1, 2, 3):
                line = fq_in.readline()
                if not line:
                    end_of_file = True
                    break
            if end_of_file:
                break

            """
            This section filters and adjusts the qualities to check. It handles cases of irregular read-lengths as well.
            """
            qualities_to_check = convert_quality_string(line.strip(), qual_offset)

            if len(qualities_to_check) != read_length:
                wrong_len += 1
                continue

            # TODO Adding this section to account for quality score "shape" in a fastq
            # shape_curves.append(qualities_to_check)

            records_read += 1

            for j in range(read_length):
                # The qualities of each read_position_scores
                quality_bin = take_closest(quality_scores, qualities_to_check[j])
                bin_index = quality_scores.index(quality_bin)
                temp_q_count[j][bin_index] += 1
                qual_score_counter[quality_bin] += 1

            if records_read % quarters == 0:
                _LOG.info(f'reading data: {(records_read / max_reads) * 100:.0f}%')

    _LOG.info(f'reading data: 100%')
    if end_of_file:
        _LOG.info(f'{records_read} records read before end of file.')
    _LOG.debug(f'{wrong_len} total reads had a length other than {read_length} ({wrong_len/max_reads:.0f}%)')

    avg_std_by_pos = []
    q_count_by_pos = np.asarray(temp_q_count)
    for i in range(read_length):
        this_counts = q_count_by_pos[i]
        expanded_counts = expand_counts(this_counts, quality_scores)
        average_q = np.average(expanded_counts)
        st_d_q = np.std(expanded_counts)
        avg_std_by_pos.append((average_q, st_d_q))

    # TODO In progress, working on ensuring the error model produces the right shape
    # shape_curves = pd.DataFrame(shape_curves)
    # columns = list(range(1, 11)) + list(range(15, len(shape_curves[0]), 5))
    # shape_curves_plot = shape_curves.iloc[columns]
    # averages = []
    # for name, value in shape_curves_plot.items():
    #     averages.append(np.average(value))
    # sns.boxplot(data=shape_curves_plot, fliersize=0)
    # plt.show()

    # Calculates the average error rate
    tot_bases = read_length * records_read
    avg_err = 0
    for score in quality_scores:
        error_val = 10. ** (-score / 10.)
        _LOG.info(f"q_score={score}, error value={error_val:e}, count={qual_score_counter[score]}")
        avg_err += error_val * (qual_score_counter[score] / tot_bases)
    _LOG.info(f'Average error rate for dataset: {avg_err}')

    # Generate the sequencing error model with default average error rate
    return avg_std_by_pos, avg_err, read_length


def plot_stuff(init_q, real_q, q_range, prob_q, actual_readlen, plot_path):
    """
    Return the factorial of n, an exact integer >= 0.

    >>> plot_stuff(init_q=30, real_q==30, q_range=0, prob_q=0.5, actual_readlen=150, plot_path=output)
    'neat/model_sequencing_error/output.svg'

    """

    plt.rcParams.update({'font.size': 14, 'font.weight': 'bold', 'lines.linewidth': 3})

    plt.figure(1)
    Z = np.array(init_q).T
    X, Y = np.meshgrid(range(0, len(Z[0]) + 1), range(0, len(Z) + 1))
    plt.pcolormesh(X, Y, Z, vmin=0., vmax=0.25)
    plt.axis([0, len(Z[0]), 0, len(Z)])
    plt.yticks(range(0, len(Z), 10), range(0, len(Z), 10))
    plt.xticks(range(0, len(Z[0]), 10), range(0, len(Z[0]), 10))
    plt.xlabel('Read Position')
    plt.ylabel('Quality Score')
    plt.title('Q-Score Prior Probabilities')
    plt.colorbar()

    plt.show()

    v_min_log = [-4, 0]
    min_val = 10 ** v_min_log[0]
    q_labels = [str(n) for n in range(q_range[0], q_range[1] + 1) if n % 5 == 0]
    print(q_labels)
    q_ticks_x = [int(n) + 0.5 for n in q_labels]
    q_ticks_y = [(real_q - int(n)) - 0.5 for n in q_labels]

    for p in range(1, actual_readlen, 10):
        current_data = np.array(prob_q[p])
        for i in range(len(current_data)):
            for j in range(len(current_data[i])):
                current_data[i][j] = max(min_val, current_data[i][j])

        # matrix indices:		pcolormesh plotting:	plot labels and axes:
        #
        #      y				   ^					   ^
        #	   -->				 x |					 y |
        #  x |					    -->					    -->
        #    v 					    y					    x
        #
        # to plot a MxN matrix 'Z' with rowNames and colNames we need to:
        #
        # pcolormesh(X,Y,Z[::-1,:])		# invert x-axis
        # # swap x/y axis parameters and labels, remember x is still inverted:
        # xlim([yMin,yMax])
        # ylim([M-xMax,M-xMin])
        # xticks()
        #

        plt.figure(p + 1)
        z = np.log10(current_data)
        x, y = np.meshgrid(range(0, len(Z[0]) + 1), range(0, len(Z) + 1))
        plt.pcolormesh(x, y, z[::-1, :], vmin=v_min_log[0], vmax=v_min_log[1], cmap='jet')
        plt.xlim([q_range[0], q_range[1] + 1])
        plt.ylim([real_q - q_range[1] - 1, real_q - q_range[0]])
        plt.yticks(q_ticks_y, q_labels)
        plt.xticks(q_ticks_x, q_labels)
        plt.xlabel('\n' + r'$Q_{i+1}$')
        plt.ylabel(r'$Q_i$')
        plt.title('Q-Score Transition Frequencies [Read Pos:' + str(p) + ']')
        cb = plt.colorbar()
        cb.set_ticks([-4, -3, -2, -1, 0])
        cb.set_ticklabels([r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$'])

    # plt.tight_layout()
    plt.show()
    plt.savefig(f'neat/model_sequencing_error/{plot_path}', format='svg')
