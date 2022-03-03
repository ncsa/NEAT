import bisect
import copy
import random
from typing import Union

import numpy as np
import pandas as pd

from source.constants_and_defaults import LOW_PROBABILITY_THRESHOLD
from source.error_handling import premature_exit, print_and_log


def mean_ind_of_weighted_list(candidate_list: list) -> int:
    """
    Returns the index of the mean of a weighted list

    :param candidate_list: weighted list
    :return: index of mean
    """
    my_mid = sum(candidate_list) / 2.0
    my_sum = 0.0
    for i in range(len(candidate_list)):
        my_sum += candidate_list[i]
        if my_sum >= my_mid:
            return i


class DiscreteDistribution:
    def __init__(self, values, weights, rng_val, degenerate_val=None):

        self.values = values
        self.weights = weights
        # some sanity checking
        if not len(weights):
            print_and_log('Probability vector given to DiscreteDistribution() are 0-length.', 'error')
            premature_exit(1)

        # if probability of all input events is 0, consider it degenerate and always return the first value
        self.degenerate = False
        if degenerate_val:
            self.degenerate = True
            self.degenerate_val = degenerate_val

    def sample(self, size=1) -> list:
        """
        This is one of the slowest parts of the code. Or it just gets hit the most times. Will need
        to investigate at some point.
        :return: Since this function is selecting an item from a list, and the list could theoretically be anything,
        then in a broad sense this function returns a list item or a generic object. But I'm fairly confident that most
        of these uses will be lists of ints or floats, but will investigate further
        """

        if self.degenerate:
            return self.degenerate_val
        elif size == 1:
            return random.choices(self.values, self.weights)[0]
        else:
            return random.choices(self.values, self.weights, k=size)


# takes k_range, lambda, [0,1,2,..], returns a DiscreteDistribution object
# with the corresponding to a poisson distribution

def poisson_list(k_range: range, input_lambda: float) -> DiscreteDistribution:
    min_weight = 1e-12
    # if input_lambda is less than some arbitrarily small number,
    # return a degenerate distribution that will always give 0.
    if input_lambda < min_weight:
        return DiscreteDistribution([1], [0], degenerate_val=0)
    # Honestly not completely sure where this code comes from or how to interpret it exactly. I think it's
    # derived from the Poisson PDF: (e^-lambda * lambda^x)/x!
    # But it's like the log of that list for x element of the k_range
    log_factorial_list = [0.0]
    for k in k_range[1:]:
        log_factorial_list.append(np.log(k) + log_factorial_list[k - 1])
    w_range = [np.exp(k * np.log(input_lambda) - input_lambda - log_factorial_list[k]) for k in k_range]
    # Here originally was the code from v2.0:
    # w_range = [n for n in w_range if n >= min_weight]
    # return DiscreteDistribution(w_range, k_range[:len(w_range)])
    # I think what was happening in that filtering step was that the Discrete Distribution was skewing lower.
    # Say input_lambda = 4.25. That means that we are shooting for roughly 4 variants. The original w_range
    # has a max value at the 4th element, and follows a distribution that makes sense, centered on 4.
    # But when the filtering happens, they then just took the first n elements of the k-range, where n = len(w_range)
    # For higher lambda values, this has the effect of skewing the weights downward. For example, an input
    # lambda of 80.7525 has a w_range.index(max(w_range)) = 80 (the index of the maximum value). This  makes sense
    # as sampling from that set should give you a number near 80. But after filtering, there are only 62 elements left
    # in the set, and the the max value moved from 80 to 28. So now when it randomly samples from this set,
    # it will pick a value near 28, which is much lower than we want.
    # My solution, then, is to used a pandas dataframe to keep the counts centered properly on the lambda value
    # during the filtering step.
    df = pd.DataFrame({'counts': k_range, 'weights': w_range})
    filtered_df = df[df.weights > min_weight].reset_index()
    if len(filtered_df) <= 1:
        return DiscreteDistribution([1], [0], degenerate_val=0)
    return DiscreteDistribution(filtered_df.weights, filtered_df.counts)


# quantize a list of values into blocks
def quantize_list(list_to_quantize):
    min_prob = 1e-12
    quant_blocks = 10
    sum_list = float(sum(list_to_quantize))
    sorted_list = sorted([n for n in list_to_quantize if n >= min_prob * sum_list])
    if len(sorted_list) == 0:
        return None
    qi = []
    for i in range(quant_blocks):
        # qi.append(sorted_list[int((i)*(len(sorted_list)/float(quant_blocks)))])
        qi.append(sorted_list[0] + (i / float(quant_blocks)) * (sorted_list[-1] - sorted_list[0]))
    qi.append(1e12)
    running_list = []
    prev_bi = None
    prev_i = None
    for i in range(len(list_to_quantize)):
        if list_to_quantize[i] >= min_prob * sum_list:
            bi = bisect.bisect(qi, list_to_quantize[i])
            # print i, l[i], qi[bi-1]
            if prev_bi is not None:
                if bi == prev_bi and prev_i == i - 1:
                    running_list[-1][1] += 1
                else:
                    running_list.append([i, i, qi[bi - 1]])
            else:
                running_list.append([i, i, qi[bi - 1]])
            prev_bi = bi
            prev_i = i
    return running_list
