import numpy as np
import random
import pickle
import copy
from source.probability import DiscreteDistribution
from source.constants_and_models import DEFAULT_MODEL_1, DEFAULT_MODEL_2, TRI_IND, \
    ALL_TRI, ALL_IND, ALLOWED_NUCL, NUC_IND
from source.error_handling import premature_exit, print_and_log


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


def parse_input_mutation_model(model=None, which_default=1):
    """
    parse mutation model pickle file
    :param model: model to read (if none, this will select a default)
    :param which_default: 1 = standard mutation model, 2 = cancer mutation model
    """

    """
    Default mutation model descriptions:
    DEFAULT_MODEL_1 = [DEFAULT_1_OVERALL_MUT_RATE = 0.001,
                       DEFAULT_1_HOMOZYGOUS_FREQ = 0.01,
                       DEFAULT_1_INDEL_FRACTION = 0.05,
                       DEFAULT_1_INS_VS_DEL = 0.6,
                       DEFAULT_1_INS_LENGTH_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                       DEFAULT_1_INS_LENGTH_WEIGHTS = [0.4, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.034, 0.033, 0.033],
                       DEFAULT_1_DEL_LENGTH_VALUES = [1, 2, 3, 4, 5],
                       DEFAULT_1_DEL_LENGTH_WEIGHTS = [0.3, 0.2, 0.2, 0.2, 0.1],
                       example_matrix_1 = [[0.0, 0.15, 0.7, 0.15],
                                           [0.15, 0.0, 0.15, 0.7],
                                           [0.7, 0.15, 0.0, 0.15],
                                           [0.15, 0.7, 0.15, 0.0]],
                       DEFAULT_1_TRI_FREQS = [copy.deepcopy(example_matrix_1) for _ in range(16)] 
                                             (i.e., 16 copies of the example matrix),
                       DEFAULT_1_TRINUC_BIAS = [1. / float(len(ALL_TRI)) for _ in ALL_TRI]
                                               (i.e., reciprocal of the number of trinucleotide combinations, 
                                               which is constant)
    """

    if which_default == 1:
        out_model = [copy.deepcopy(n) for n in DEFAULT_MODEL_1]
    elif which_default == 2:
        out_model = [copy.deepcopy(n) for n in DEFAULT_MODEL_2]
    else:
        print_and_log('Unknown default mutation model specified.', 'error')
        premature_exit(1)

    if model is not None:
        pickle_dict = pickle.load(open(model, "rb"))
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

        self.debug = options.debug

        # Lead mutation models
        self.mutation_model = parse_input_mutation_model(options.mutation_model)
        if options.cancer:
            self.cancer_model = parse_input_mutation_model(options.cancer_model, 2)

        if self.debug:
            print_and_log("Mutation models loaded", 'debug')

        # Load sequencing error model



