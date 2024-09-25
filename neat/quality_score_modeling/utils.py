import pysam
import pandas as pd
import numpy as np
from scipy.stats import norm
import pathlib
import pickle
from Bio import SeqIO

__all__ = [
    "make_qual_score_list",
    "apply_markov_chain",
    "save_file"
]

def make_qual_score_list(fastq_file=None, bam_file=None):
    """Takes an input FASTQ or BAM file and creates lists of quality scores. This becomes a DataFrame, which
    will be pre-processed for Markov chain analysis.

    If both FASTQ and BAM are provided, BAM file will be used.

    :param fastq_file: Path to the input FASTQ file (optional if BAM is provided)
    :param bam_file: Path to the input BAM file (optional if FASTQ is provided)
    :return: DataFrame with quality scores.
    """

    if bam_file:

        index = f"{bam_file}.bai"

        if not pathlib.Path(index).exists():
            print("No index found, creating one.")

            pysam.index(bam_file)

        file_to_parse = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
        num_recs = file_to_parse.count()
        print(f"{num_recs} records to parse")

        modulo = round(num_recs / 9)

        qual_list = []
        i = 0
        j = 0

        def print_update(number, factor, percent):

            if number % factor == 0:
                percent += 10
                print(f"{percent}% complete", end="\r")

            return percent

        print("Parsing file")

        for item in file_to_parse.fetch():

            if item.is_unmapped or "S" in item.cigarstring:
                i += 1
                j = print_update(i, modulo, j)

                continue

            # mapping quality scores

            align_qual = item.query_alignment_qualities

            # append to master lists

            qual_list.append(align_qual)
            i += 1
            j = print_update(i, modulo, j)

        print(f"100% complete")
        file_to_parse.close()

    elif fastq_file:
        # Process FASTQ file if provided
        print("Parsing FASTQ file")

        qual_list = []
        for record in SeqIO.parse(fastq_file, "fastq"):
            qual_list.append(record.letter_annotations["phred_quality"])

    else:
        raise ValueError("Either fastq_file or bam_file must be provided.")

    quality_df = pd.DataFrame(qual_list)  # turn list of lists into a dataframe
    quality_df = quality_df.fillna(0)  # pre-processing - fill in missing data

    # filters to process outliers

    quality_df[quality_df > 40] = 40
    quality_df[quality_df < 0] = 0

    return quality_df


def estimate_transition_probabilities(std_dev):
    """Takes a standard deviation as an input and generates the transition probabilities with a normal
    distribution that can be used to represent a Markov process."""

    # define the probabilities for transition states based on a normal distribution

    transition_probs = {
        -3: norm.pdf(-3, 0, std_dev),
        -2: norm.pdf(-2, 0, std_dev),
        -1: norm.pdf(-1, 0, std_dev),
        0: norm.pdf(0, 0, std_dev),
        1: norm.pdf(1, 0, std_dev),
        2: norm.pdf(2, 0, std_dev),
        3: norm.pdf(3, 0, std_dev),
    }

    # normalize the probabilities to sum to 1

    total_prob = sum(transition_probs.values())

    for k in transition_probs:
        transition_probs[k] /= total_prob

    return transition_probs


def apply_markov_chain(quality_df, noise_level=10, std_dev=2):
    """Takes a data frame representing quality scores by position along a read and parameters to increase
    variability in the Markov process as inputs and generates predictions based on an ergodic Markov chain.
    Generates a data frame with simulated reads."""

    transition_probs = estimate_transition_probabilities(std_dev)
    num_rows, num_cols = quality_df.shape

    count = 0
    markov_preds = []

    for row in quality_df.iterrows():

        qualities = row[1].values
        pred_qualities = np.zeros_like(qualities)
        pred_qualities[0] = qualities[0]  # initial state

        # print(count, ":", qualities)
        row_mean = np.mean(qualities)
        row_median = np.median(qualities)
        row_std = np.std(qualities)

        for i in range(1, len(quality_df.columns)):
            prev_quality = pred_qualities[i - 1]
            transitions = list(transition_probs.keys())
            probabilities = list(transition_probs.values())
            next_quality = np.random.choice(transitions, p=probabilities)

            pred_qualities[i] = max(0, prev_quality + next_quality)  # ensuring no negative qualities

            # the noise parameter prevents long stretches of the predicted quality scores being very similar

            pred_qualities[i] += np.random.normal(0, noise_level)  # add some noise to the predictions
            pred_qualities[i] = min(max(pred_qualities[i], 0), 40)  # finalize range

        # print(count, "mean:", row_mean, "median:", row_median, "st dev:", row_std)
        count += 1

        for i in range(1, len(quality_df.columns)):

            if pred_qualities[i] < row_mean - 2 * row_std:

                if np.random.rand() < 0.95:  # 95% chance to substitute abnormally low quality scores

                    # uses median and standard deviation from read (not the mean because of outliers)

                    new_quality = np.random.normal(row_median, row_std)
                    pred_qualities[i] = min(max(new_quality, 0), 40)

        # the maximum predicted quality score should be derived from the data

        max_quality = np.max(qualities)
        pred_qualities = np.clip(pred_qualities, 0, max_quality)

        markov_preds.append(pred_qualities)

        # randomly sample 30% of the total quality scores for a given read to have the maximum value

        num_samples = int(0.3 * len(quality_df.columns))  # TO DO: make this a parameter
        sample_indices = np.random.choice(len(quality_df.columns), num_samples, replace=False)
        pred_qualities[sample_indices] = max_quality

    markov_preds_df = pd.DataFrame(markov_preds)

    # apply final linear transformations

    edge_len = int(len(quality_df.columns) * 0.01)
    # mid_start = int(len(quality_df.columns) * 0.40)
    # mid_end = int(len(quality_df.columns) * 0.60)

    markov_preds_df.iloc[:, :edge_len] -= 5
    markov_preds_df.iloc[:, :edge_len] = markov_preds_df.iloc[:, :edge_len].clip(lower=0)

    markov_preds_df.iloc[:, -edge_len:] -= 5
    markov_preds_df.iloc[:, -edge_len:] = markov_preds_df.iloc[:, -edge_len:].clip(lower=0)

    # markov_preds_df.iloc[:, mid_start:mid_end] += 2

    return markov_preds_df


def save_file(df, csv_file_path, pickle_file_path):
    """Saves the dataframe to a CSV file and a pickle file."""

    # df.to_csv(csv_file_path, index=False) # figure out how to turn off by default

    with open(pickle_file_path, "wb") as f:
        pickle.dump(df, f)

    print(f"Data saved to {csv_file_path} and {pickle_file_path}")


# def load_markov_predictions(pickle_file):
#     with open(pickle_file, "rb") as f:
#         markov_preds_df = pickle.load(f)
#     return markov_preds_df
