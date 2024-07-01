import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
from scipy.stats import ttest_ind, ttest_rel, f_oneway, norm, levene, shapiro
from sklearn.utils import resample
import pathlib
import pickle


def make_qual_score_list(bam_file):
    """Takes an input BAM file and creates lists of quality scores. This becomes a data frame, which will be
    pre-processed for Markov chain analysis."""

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

        print(count, ":", qualities)
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

        print(count, "mean:", row_mean, "median:", row_median, "st dev:", row_std)
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


def plot_heatmap(y_preds_df, file_path):
    """Takes a dataframe of predicted quality scores and plots a seaborn heatmap to visualize them."""

    sns.heatmap(y_preds_df, vmin=0, vmax=max(markov_preds_df.max()), cmap="viridis")
    plt.savefig(file_path)
    print("Heatmap plotted")


def save_file(df, csv_file_path, pickle_file_path):
    """Saves the dataframe to a CSV file and a pickle file."""

    df.to_csv(csv_file_path, index=False)

    with open(pickle_file_path, "wb") as f:
        pickle.dump(df, f)

    print(f"Data saved to {csv_file_path} and {pickle_file_path}")


def compare_quality_scores(quality_df, markov_preds_df):
    """Compares the means and variances of quality scores within each row of the original and predicted data frames."""

    mean_p_values = []
    variance_p_values = []

    for i in range(len(quality_df)):
        original = quality_df.iloc[i].values
        predicted = markov_preds_df.iloc[i].values

        original_mean = np.mean(original)
        predicted_mean = np.mean(predicted)

        # print the pairs of means

        print(f'Row {i}: Original Mean = {original_mean}, Predicted Mean = {predicted_mean}')

        # test for equality of means

        t_stat, mean_p_value = ttest_ind(original, predicted, equal_var=False)
        mean_p_values.append(mean_p_value)

        # test for equality of variances

        stat, variance_p_value = levene(original, predicted)
        variance_p_values.append(variance_p_value)

    return mean_p_values, variance_p_values


def plot_comparison_results(mean_p_values, variance_p_values):
    """Plots the comparison results for means and variances of quality scores."""

    rows = range(len(mean_p_values))

    fig, ax = plt.subplots(1, 2, figsize=(12, 6))

    ax[0].scatter(rows, mean_p_values, marker='o')
    ax[0].set_title('P-values for Equality of Means')
    ax[0].set_xlabel('Row Index')
    ax[0].set_ylabel('P-value')

    ax[1].scatter(rows, variance_p_values, marker='o')
    ax[1].set_title('P-values for Equality of Variances')
    ax[1].set_xlabel('Row Index')
    ax[1].set_ylabel('P-value')

    plt.tight_layout()
    plt.show()


def calculate_row_means(df):
    return df.mean(axis=1)


def calculate_row_variances(df):
    return df.var(axis=1)


def bootstrap_p_values(list1, list2, n_bootstrap=1000):
    bootstrapped_p_values = []

    for _ in range(n_bootstrap):
        sample1 = resample(list1, replace=True)
        sample2 = resample(list2, replace=True)

        t_stat, p_value = ttest_ind(sample1, sample2, equal_var=False)
        bootstrapped_p_values.append(p_value)

    return bootstrapped_p_values


def permutation_test(list1, list2, n_permutations=10000):
    observed_diff = abs(np.mean(list1) - np.mean(list2))
    combined = np.concatenate([list1, list2])
    perm_diffs = np.zeros(n_permutations)

    for i in range(n_permutations):
        np.random.shuffle(combined)
        perm_list1 = combined[:len(list1)]
        perm_list2 = combined[len(list1):]
        perm_diffs[i] = abs(np.mean(perm_list1) - np.mean(perm_list2))

    p_value = np.sum(perm_diffs >= observed_diff) / n_permutations
    return p_value


def test_normality(quality_df, markov_preds_df):
    """Tests normality of quality scores within each row using the Shapiro-Wilk test."""

    quality_df_normality = []
    markov_preds_df_normality = []

    for i in range(len(quality_df)):
        original = quality_df.iloc[i].values
        predicted = markov_preds_df.iloc[i].values

        # test normality for the original data

        stat, p_value = shapiro(original)
        quality_df_normality.append(p_value > 0.05)

        # test normality for the predicted data

        stat, p_value = shapiro(predicted)
        markov_preds_df_normality.append(p_value > 0.05)

    return quality_df_normality, markov_preds_df_normality


def plot_normality_results(quality_df_normality, markov_preds_df_normality):
    """Plots the normality results for quality scores in side by side heatmaps."""

    fig, ax = plt.subplots(1, 2, figsize=(12, 6))

    quality_df_normality = np.array(quality_df_normality).reshape(-1, 1)
    markov_preds_df_normality = np.array(markov_preds_df_normality).reshape(-1, 1)

    # create a color map for the heatmap

    cmap = sns.color_palette(["orange", "blue"])

    sns.heatmap(quality_df_normality, ax=ax[0], cbar=False, cmap=cmap)
    ax[0].set_title('Normality of Original Quality Scores')
    ax[0].set_xlabel('Quality Scores')
    ax[0].set_ylabel('Row Index')

    sns.heatmap(markov_preds_df_normality, ax=ax[1], cbar=False, cmap=cmap)
    ax[1].set_title('Normality of Predicted Quality Scores')
    ax[1].set_xlabel('Quality Scores')

    legend_elements = [Patch(facecolor='orange', edgecolor='black', label='Non-normal'),
                       Patch(facecolor='blue', edgecolor='black', label='Normal')]
    fig.legend(handles=legend_elements, loc='lower right', title='Distribution')

    plt.tight_layout()
    plt.show()


# example usage

# bam_file = "/Users/keshavgandhi/Downloads/H1N1.bam"

bam_file = "/Users/keshavgandhi/Downloads/subsample_3.125.bam"
quality_df = make_qual_score_list(bam_file)

markov_preds_df = apply_markov_chain(quality_df)

# plot_heatmap(markov_preds_df, 'markov_chain_heatmap.svg')
# save_to_csv_and_pickle(markov_preds_df, 'markov_preds.csv', 'markov_preds.pickle')

sns.heatmap(quality_df, vmin=0, vmax=max(quality_df.max()), cmap='viridis')
sns.heatmap(markov_preds_df, vmin=0, vmax=max(markov_preds_df.max()), cmap='viridis')

# markov_preds_df

# for i in range (1, max(markov_preds_df)):
#     print(max(markov_preds_df[i]))

# quality_df.iloc[100][25:75]
# markov_preds_df.iloc[100][25:75]

bam_file = "/Users/keshavgandhi/Downloads/H1N1.bam"
test_df = make_qual_score_list(bam_file)
markov_preds_df = apply_markov_chain(test_df)

# compare quality scores

mean_p_values, variance_p_values = compare_quality_scores(test_df, markov_preds_df)

markov_means = calculate_row_means(markov_preds_df).tolist()
quality_means = calculate_row_means(test_df).tolist()

markov_variances = calculate_row_variances(markov_preds_df).tolist()
quality_variances = calculate_row_variances(test_df).tolist()

# perform permutation test

permutation_p_value_means = permutation_test(markov_means, quality_means)

# perform two-sample t-test

t_stat_means, ttest_p_value_means = ttest_ind(markov_means, quality_means, equal_var=False)

# bootstrap analysis for means

bootstrapped_p_values_means = bootstrap_p_values(markov_means, quality_means)
mean_bootstrap_p_value_means = np.mean(bootstrapped_p_values_means)
std_bootstrap_p_value_means = np.std(bootstrapped_p_values_means)

print(f'Permutation test p-value (means): {permutation_p_value_means}')
print(f'Two-sample t-test p-value (means): {ttest_p_value_means}')
print(f'Bootstrap mean p-value (means): {mean_bootstrap_p_value_means}')
print(f'Bootstrap p-value standard deviation (means): {std_bootstrap_p_value_means}')

# perform permutation test for variances

permutation_p_value_variances = permutation_test(markov_variances, quality_variances)

# perform two-sample t-test for variances

t_stat_variances, ttest_p_value_variances = ttest_ind(markov_variances, quality_variances, equal_var=False)

# perform bootstrap analysis for variances

bootstrapped_p_values_variances = bootstrap_p_values(markov_variances, quality_variances)
mean_bootstrap_p_value_variances = np.mean(bootstrapped_p_values_variances)
std_bootstrap_p_value_variances = np.std(bootstrapped_p_values_variances)

print(f'Permutation test p-value (variances): {permutation_p_value_variances}')
print(f'Two-sample t-test p-value (variances): {ttest_p_value_variances}')
print(f'Bootstrap mean p-value (variances): {mean_bootstrap_p_value_variances}')
print(f'Bootstrap p-value standard deviation (variances): {std_bootstrap_p_value_variances}')

# plot comparison results

plot_comparison_results(mean_p_values, variance_p_values)

# test normality

quality_df_normality_test, markov_preds_df_normality_test = test_normality(test_df, markov_preds_df)

# plot normality results

plot_normality_results(quality_df_normality_test, markov_preds_df_normality_test)
