import numpy as np
import random
import pickle
import copy
import gzip
from source.probability import DiscreteDistribution
from source.constants_and_models import DEFAULT_MUTATION_MODEL, DEFAULT_CANCER_MUTATION_MODEL, TRI_IND, \
    ALL_TRI, ALL_IND, ALLOWED_NUCL, NUC_IND
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
    see constants_and_models.py for a complete description of the default models.
    :param model: model to read (if none, this will select a default)
    :param is_cancer: False = standard mutation model, True = cancer mutation model
    """

    if not is_cancer:
        out_model = [copy.deepcopy(n) for n in DEFAULT_MUTATION_MODEL]
    else:
        out_model = [copy.deepcopy(n) for n in DEFAULT_CANCER_MUTATION_MODEL]

    if model:
        pickle_dict = pickle.load(gzip.open(model, "rb"))
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

        # Load sequencing error model
        mssg = "Problem loading Sequencing error model data @error_model"
        self.sequencing_error_model = pickle_load_model(gzip.open(options.error_model, 'rb'), mssg)

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
