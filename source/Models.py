import copy
import gzip
import pickle
import math

import numpy as np
import pandas as pd

from source.constants_and_defaults import DEFAULT_MUTATION_MODEL, DEFAULT_CANCER_MUTATION_MODEL, TRI_IND, \
    ALL_TRI, ALL_IND, NUC_IND, ALLOWED_NUCL
from source.error_handling import premature_exit, print_and_log
from source.probability import DiscreteDistribution, mean_ind_of_weighted_list


def pickle_load_model(file, mssg) -> list:
    """
    This will pickle load a file for processing, checking some standard errors that can happen with the pickle
    files created by the utilities.
    :file: the pathlib object to the file.
    :mssg: the error message to print if things go awry
    :return: returns the contents if the file, which should all be lists
    """
    try:
        return pickle.load(file)
    except IOError as e:
        print_and_log(mssg, 'error')
        print_and_log(str(e), 'error')
        premature_exit(1)
    except EOFError as e:
        print_and_log(mssg, 'error')
        print_and_log(str(e), 'error')
        premature_exit(1)
    except ValueError as e:
        print_and_log(mssg, 'error')
        print_and_log(str(e), 'error')
        premature_exit(1)


def parse_input_mutation_model(model, is_cancer: bool = False):
    """
    parse mutation model pickle file
    see constants_and_defaults.py for a complete description of the default models.
    :param model: model to read (if none, this will select a default)
    :param is_cancer: False = standard mutation model, True = cancer mutation model
    """

    if not is_cancer:
        out_model = [copy.deepcopy(n) for n in DEFAULT_MUTATION_MODEL]
    else:
        out_model = [copy.deepcopy(n) for n in DEFAULT_CANCER_MUTATION_MODEL]

    if model:
        mssg = "Problem loading the mutation model."
        pickle_dict = pickle_load_model(gzip.open(model, "rb"), mssg)
        out_model[0] = pickle_dict['AVG_MUT_RATE']
        out_model[2] = 1. - pickle_dict['SNP_FREQ']

        ins_list = pickle_dict['INDEL_FREQ']
        if len(ins_list):
            ins_count = sum([ins_list[k] for k in ins_list.keys() if k >= 1])
            del_count = sum([ins_list[k] for k in ins_list.keys() if k <= -1])
            ins_vals = [k for k in sorted(ins_list.keys()) if k >= 1]
            ins_weight = [ins_list[k] / float(ins_count) for k in ins_vals]
            del_vals = [k for k in sorted([abs(k) for k in ins_list.keys() if k <= -1])]
            del_weight = [ins_list[-k] / float(del_count) for k in del_vals]
        else:  # degenerate case where no indel stats are provided
            ins_count = 1
            del_count = 1
            ins_vals = [1]
            ins_weight = [1.0]
            del_vals = [1]
            del_weight = [1.0]
        out_model[3] = ins_count / float(ins_count + del_count)
        out_model[4] = ins_vals
        out_model[5] = ins_weight
        out_model[6] = del_vals
        out_model[7] = del_weight

        trinuc_trans_prob = pickle_dict['TRINUC_TRANS_PROBS']
        for k in sorted(trinuc_trans_prob.keys()):
            my_ind = TRI_IND[k[0][0] + k[0][2]]
            (k1, k2) = (NUC_IND[k[0][1]], NUC_IND[k[1][1]])
            out_model[8][my_ind][k1][k2] = trinuc_trans_prob[k]
        for i in range(len(out_model[8])):
            for j in range(len(out_model[8][i])):
                for l in range(len(out_model[8][i][j])):
                    # if trinuc not present in input mutation model, assign it uniform probability
                    if float(sum(out_model[8][i][j])) < 1e-12:
                        out_model[8][i][j] = [0.25, 0.25, 0.25, 0.25]
                    else:
                        out_model[8][i][j][l] /= float(sum(out_model[8][i][j]))

        trinuc_mut_prob = pickle_dict['TRINUC_MUT_PROB']
        which_have_we_seen = {n: False for n in ALL_TRI}
        trinuc_mean = np.mean(list(trinuc_mut_prob.values()))
        for trinuc in trinuc_mut_prob.keys():
            out_model[9][ALL_IND[trinuc]] = trinuc_mut_prob[trinuc]
            which_have_we_seen[trinuc] = True
        for trinuc in which_have_we_seen.keys():
            if not which_have_we_seen[trinuc]:
                out_model[9][ALL_IND[trinuc]] = trinuc_mean

    return out_model


def score_breakdown(quality_bins: list):
    """
    breaks down the quality scores into a dicitonary of the interval and error rate for that bin.
    This code requires quality bins to be a monotonic increasing list of scores corresponding to the phred score.
    """
    # This block of code breaks down the quality scores into intervals for analysis.
    score_breakdown = {}
    low = 0
    mid_point = (quality_bins[1] - quality_bins[0]) // 2
    high = quality_bins[0] + mid_point
    # This is the inverse of the Phred quality score formula f(x) = -10 * log10(x)
    error_rate = 10. ** (-quality_bins[0] / 10.)
    # We'll precompute, for each score, both the interval it spans and the error_rate it represents
    score_breakdown[quality_bins[0]] = {'interval': (low, high), 'error_rate': error_rate}
    low = high + 1

    for i in range(1, len(quality_bins) - 1):
        mid_point = (quality_bins[i + 1] - quality_bins[i]) // 2
        high = quality_bins[i] + mid_point
        error_rate = 10. ** (-quality_bins[i] / 10.)
        score_breakdown[quality_bins[i]] = {'interval': (low, high), 'error_rate': error_rate}
        low = high + 1

    error_rate = 10. ** (-quality_bins[-1] / 10.)
    score_breakdown[quality_bins[-1]] = {'interval': (low, math.inf), 'error_rate': error_rate}
    return score_breakdown


class QualityScoreModel:
    """
    This class represents the parameters needed to model quality scores.
    """
    def __init__(self, options):
        """
        This is a function version of the old SequencingErrors class. May convert to a class
        if it seems valuable. For now, let's try this.

        New input model:

        Dictionary of the probabilities and values:
            - quality_scores: possible values of quality scores. This is a dictionary of the form:
                - The keys are the possible bins for the quality scores. This is input when the model is being made
                - The values are themselves dictionaries:
                    - the keys are "interval" and "error_rate"
                    - the value for "interval" is a tuple, which gives the high and low value of the range the
                      bin represents
                    - the value for "error_rate" is just a float that is the inverse of the phred score the bin
                      represents.
            - quality_score_probabilities: pandas dataframe indexed by position along the read, with each row
                                           representing the probabilities of each quality bin at that location.
            - quality_offset: a number representing the offset used to translate from quality score to ascii character,
                              usually this is 33, but some machines differ.
            - average_error: the estimated average error of the entire dataset
            - is_uniform: a boolean that indicates if the quality scores have uniform error.
            - error_parameters:
                - error_params[0] - A nucleotide substitution matrix
                - error_params[1] - The chance that he error is an indel
                - error_params[2] - The probability distribution of the possible lengths of indels
                - error_params[3] - The possible lengths of indels
                - error_params[4] - The probability of an indel being an insertion
                - error_params[5] - Probability distribution for the 4 nucleotides

        """
        mssg = "Problem loading Sequencing error model data @error_model"
        error_model = pickle_load_model(gzip.open(options.error_model, 'rb'), mssg)

        self.quality_scores = error_model['quality_scores'].keys()
        self.quality_score_literal = {x: y['error_rate'] for x, y in error_model['quality_scores'].items()}
        self.intervals = {x: y['interval'] for x, y in error_model['quality_scores'].items()}
        self.average_error = error_model['average_error']
        if error_model['is_uniform']:
            probs = pd.DataFrame({row: self.average_error for x in error_model['quality_scores']}, columns=error_model['quality_scores'])
        else:
            probs = error_model['quality_score_probabilities']
        self.quality_score_probabilities = probs
        self.substitution_matrix = error_model['error_parameters'][0]
        # Probability that if there is an error, it is an indel (v snp)
        self.indel_probability = error_model['error_parameters'][1]
        # Possible lengths for indels, along with the weights of each
        self.indel_length_model = dict(zip(error_model['error_parameters'][3], error_model['error_parameters'][2]))
        # Nucleotides and their probability values for a random choice
        self.nucleotide_probabilities = dict(zip(ALLOWED_NUCL, error_model['error_parameters'][5]))
        self.quality_offset = error_model['quality_offset']

        # Not sure what effect rescaling the quality scores will have on binned qualities.
        self.error_scale = 1.0
        if options.rescale_qualities:
            self.error_scale = options.avg_seq_error / error_model['average_error']

        # The point of this is supposed to be spread out the error model across a broader read length,
        # if the read length the user requests is different from what the sequencing error model calculated.
        self.quality_index_remap = range(options.read_len)
        if options.read_len != len(error_model['quality_score_probabilities']):
            print_and_log(f'Read length of error model ({error_model["quality_score_probabilities"]}) '
                          f'does not match -R value ({options.read_len}), rescaling model...', 'warning')
            self.quality_index_remap = [max([1, len(error_model["quality_score_probabilities"]) * n //
                                             options.read_len]) for n in range(options.read_len)]

    def get_sequencing_errors(self, read_data, is_reverse_strand=False):
        """
        Inserts errors of type substitution, insertion, or deletion into read_data, and assigns a quality score
        based on the container model.
        :param read_data: sequence to insert errors into
        :param is_reverse_strand: whether to treat this as the reverse strand or not
        :return: modified sequence and associated quality scores
        """



class Models:
    """
    Sets up models for use in the rest of the program.
    """
    def __init__(self, options):

        # Lead mutation models
        self.mutation_model = parse_input_mutation_model(options.mutation_model)

        self.cancer_model = None
        if options.cancer:
            self.cancer_model = parse_input_mutation_model(options.cancer_model, True)

        if options.debug:
            print_and_log("Mutation models loaded", 'debug')

        # parse sequencing error file
        self.quality_score_model = QualityScoreModel(options)

        if options.debug:
            print_and_log('Sequencing error model loaded', 'debug')

        mssg = "f'ERROR: problem reading @gc_model. Please check file path and try again. " \
               "This file should be the output of compute_gc.py'"
        self.gc_model = pickle_load_model(gzip.open(options.gc_model, 'rb'), mssg)

        if options.debug:
            print_and_log('GC Bias model loaded', 'debug')

        self.fraglen_model = None

        if options.paired_ended:
            if options.fragment_model:
                mssg = 'Problem loading the empirical fragment length model @fragment_model. Please check file and try' \
                       'again.'
                print_and_log("Using empirical fragment length distribution", 'info')
                potential_values, potential_prob = pickle_load_model(open(options.fragment_model, 'rb'), mssg)

                fraglen_values = []
                fraglen_probability = []
                for i in range(len(potential_values)):
                    if potential_values[1] > options.read_len:
                        fraglen_values.append(potential_values[i])
                        fraglen_probability.append(potential_prob[i])

                self.fraglen_model = DiscreteDistribution(fraglen_probability, fraglen_values)
                options.set_value('fragment_mean', fraglen_values[mean_ind_of_weighted_list(fraglen_probability)])

            # Using artificial fragment length distribution, if the parameters were specified
            # fragment length distribution: normal distribution that goes out to +- 6 standard deviations
            else:
                print_and_log(f'Using artificial fragment length distribution.', 'info')
                if options.fragment_st_dev == 0:
                    self.fraglen_model = DiscreteDistribution([1], [options.fragment_mean],
                                                              degenerate_val=options.fragment_mean)
                else:
                    potential_values = range(max(0, int(options.fragment_mean - 6 * options.fragment_st_dev)),
                                             int(options.fragment_mean + 6 * options.fragment_st_dev) + 1)
                    fraglen_values = []
                    for i in range(len(potential_values)):
                        if potential_values[i] > options.read_len:
                            fraglen_values.append(potential_values[i])
                    fraglen_probability = [np.exp(-(((n - options.fragment_mean) ** 2) /
                                                    (2 * (options.fragment_st_dev ** 2)))) for n in
                                           fraglen_values]
                    self.fraglen_model = DiscreteDistribution(fraglen_probability, fraglen_values)
            if options.debug:
                print_and_log(f'Loaded paired-end models', 'debug')