import pysam
import numpy as np
import pandas as pd
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm

def make_qual_score_list(bam_file):
    '''Takes an input BAM file and creates lists of quality scores. This becomes a data frame, which will be
    pre-processed for Markov chain analysis.'''

    index = f'{bam_file}.bai'

    if not pathlib.Path(index).exists():
        print('No index found, creating one.')

        pysam.index(bam_file)

    file_to_parse = pysam.AlignmentFile(bam_file, 'rb', check_sq=False)
    num_recs = file_to_parse.count()
    print(f'{num_recs} records to parse')

    modulo = round(num_recs / 9)

    qual_list = []
    i = 0
    j = 0

    def print_update(number, factor, percent):
        if number % factor == 0:
            percent += 10
            print(f'{percent}% complete', end='\r')
        return percent

    print('Parsing file')

    for item in file_to_parse.fetch():
        if item.is_unmapped or len(item.seq) != 249 or 'S' in item.cigarstring:
            i += 1
            j = print_update(i, modulo, j)
            continue

        # Mapping quality scores

        align_qual = item.query_alignment_qualities

        # Append to master lists

        qual_list.append(align_qual)
        i += 1
        j = print_update(i, modulo, j)

    print(f'100% complete')
    file_to_parse.close()

    # Turn list of lists into a dataframe

    quality_df = pd.DataFrame(qual_list)

    # Pre-processing - fill in missing data and sample 1% of reads

    quality_df = quality_df.fillna(0)
    quality_df = quality_df.sample(frac=0.01, axis=0, random_state=42)

    return quality_df

def estimate_transition_probabilities():

    # Define the probabilities for transition states based on a normal distribution

    std_dev = 1
    transition_probs = {
        -3: norm.pdf(-3, 0, std_dev),
        -2: norm.pdf(-2, 0, std_dev),
        -1: norm.pdf(-1, 0, std_dev),
        0: norm.pdf(0, 0, std_dev),
        1: norm.pdf(1, 0, std_dev),
        2: norm.pdf(2, 0, std_dev),
        3: norm.pdf(3, 0, std_dev)
    }

    # Normalize the probabilities to sum to 1

    total_prob = sum(transition_probs.values())
    for k in transition_probs:
        transition_probs[k] /= total_prob

    return transition_probs

def apply_markov_chain(quality_df, L=249):
    
    transition_probs = estimate_transition_probabilities()
    num_rows, num_cols = quality_df.shape

    markov_preds = []

    for row in quality_df.iterrows():
        qualities = row[1].values
        pred_qualities = np.zeros_like(qualities)
        pred_qualities[0] = qualities[0] # initial state
 
        for i in range(1, L):
            prev_quality = pred_qualities[i - 1]
            transitions = list(transition_probs.keys())
            probabilities = list(transition_probs.values())
            next_quality = np.random.choice(transitions, p=probabilities)
            pred_qualities[i] = max(0, prev_quality + next_quality) # ensuring no negative qualities

        markov_preds.append(pred_qualities)

    markov_preds_df = pd.DataFrame(markov_preds)

    # Apply final transformations

    edge_len = int(L * 0.05)
    mid_start = int(L * 0.40)
    mid_end = int(L * 0.60)

    markov_preds_df.iloc[:, :edge_len] -= 5
    markov_preds_df.iloc[:, :edge_len] = markov_preds_df.iloc[:, :edge_len].clip(lower=0)

    markov_preds_df.iloc[:, -edge_len:] -= 5
    markov_preds_df.iloc[:, -edge_len:] = markov_preds_df.iloc[:, -edge_len:].clip(lower=0)

    markov_preds_df.iloc[:, mid_start:mid_end] += 1

    return markov_preds_df

def plot_heatmap(y_preds_df, file_path):
    '''Takes a dataframe of predicted quality scores and plots a seaborn heatmap to visualize them.'''

    sns.heatmap(y_preds_df, vmin=0, vmax=40, cmap='viridis')
    plt.savefig(file_path)
    print('Heatmap plotted')

# Example usage

bam_file = '/projects/bclt/neat_data/H1N1_new.bam'
test_df = make_qual_score_list(bam_file)
markov_preds_df = apply_markov_chain(test_df)
plot_heatmap(markov_preds_df, 'markov_chain_heatmap.svg')

