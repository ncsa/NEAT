import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
from scipy.stats import ttest_ind, levene, shapiro
from sklearn.utils import resample


def plot_heatmap(y_preds_df, file_path):
    """Takes a dataframe of predicted quality scores and plots a seaborn heatmap to visualize them."""

    sns.heatmap(y_preds_df, vmin=0, vmax=max(y_preds_df.max()), cmap="viridis")
    plt.savefig(file_path)
    print("Heatmap plotted")


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
