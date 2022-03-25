import copy
import gzip
import pickle
import math
import random

import numpy as np
import pandas as pd

from bisect import bisect_left

from source.constants_and_defaults import DEFAULT_MUTATION_MODEL, DEFAULT_CANCER_MUTATION_MODEL, TRI_IND, \
    ALL_TRI, ALL_IND, NUC_IND, ALLOWED_NUCL
from source.error_handling import premature_exit, log_mssg
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
        log_mssg(mssg, 'error')
        log_mssg(str(e), 'error')
        premature_exit(1)
    except EOFError as e:
        log_mssg(mssg, 'error')
        log_mssg(str(e), 'error')
        premature_exit(1)
    except ValueError as e:
        log_mssg(mssg, 'error')
        log_mssg(str(e), 'error')
        premature_exit(1)


def parse_input_mutation_model(model, is_cancer: bool = False):
    """
    parse mutation model pickle file
    see constants_and_defaults.py for a complete description of the default models.
    :param model: model to read (if none, this will select a default)
    :param is_cancer: False = standard mutation model, True = cancer mutation model
    :return: The return is a dictionary of values for sampling mutations during the run
             avg_mut_rate = average mutation rate for the modeled data set
             homozygous_frequency = How frequently reads are homozygous (hard coded as 0.010)
                                    (note, we need to investigate the equivalent with polyploid organisms)
             indel_freq = probability that if there is a mutation, it is and indel (versus snp)
             indel_insert_percentage = probability that an indel is an insert (versus deletion)
             insert_length_values = Potential lengths of inserts
             insert_length_weights = Weights of the list above
             deletion_length_values = Potential lengths of deletions
             deletion_length_weights = Weights of the list above
             trinuc_freqs = The calculated substiution matrix for each trinucleotide
             trinuc_bias = The bias for each trinuc, relative to the others
    """

    if not is_cancer:
        # as the mutation rate evolves, this can be expanded or changed
        keys = ['avg_mut_rate',
                'homozygous_freq',
                'indel_freq',
                'indel_insert_percentage',
                'insert_length_values',
                'insert_length_weights',
                'deletion_length_values',
                'deletion_length_weights',
                'trinuc_freqs',
                'trinuc_bias']

        out_model = dict(zip(keys, DEFAULT_MUTATION_MODEL))
    else:
        keys = ['avg_mut_rate',
                'homozygous_freq',
                'indel_freq',
                'indel_insert_percentage',
                'insert_length_values',
                'insert_length_weights',
                'deletion_length_values',
                'deletion_length_weights',
                'trinuc_freqs',
                'trinuc_bias']
        out_model = dict(zip(keys, DEFAULT_CANCER_MUTATION_MODEL))

    if model:
        mssg = "Problem loading the mutation model."
        pickle_dict = pickle_load_model(gzip.open(model, "rb"), mssg)
        out_model['avg_mut_rate'] = pickle_dict['AVG_MUT_RATE']
        out_model['indel_freq'] = 1. - pickle_dict['SNP_FREQ']

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
        out_model['indel_insert_percentage'] = ins_count / float(ins_count + del_count)
        out_model['insert_length_model'] = DiscreteDistribution(ins_vals, np.array(ins_weight))
        out_model['deletion_length_model'] = DiscreteDistribution(del_vals, np.array(del_weight))

        out_model['trinuc_bias'] = pickle_dict['TRINUC_MUT_PROB']
        # Order is ACTG. Start by assuming they're all uniform
        out_model['trinuc_trans_prob'] = {x: {'A': 0.25, 'C': 0.25, 'T': 0.25, 'G': 0.25} for x in ALL_TRI}
        for transition in pickle_dict['TRINUC_TRANS_PROBS']:
            goes_to = transition[1][1]
            out_model['trinuc_trans_prob'][transition[0]][goes_to] = pickle_dict['TRINUC_TRANS_PROBS'][transition]

    else:
        # These will just be uniform
        out_model['trinuc_mut_prob'] = DiscreteDistribution(ALL_TRI, [1] * len(ALL_TRI))
        out_model['trinuc_trans_prob'] = {x: {'A': 0.25, 'C': 0.25, 'T': 0.25, 'G': 0.25} for x in ALL_TRI}

        out_model['insert_length_model'] = DiscreteDistribution(out_model['insert_length_values'],
                                                                out_model['insert_length_weights'])
        out_model['deletion_length_model'] = DiscreteDistribution(out_model['deletion_length_values'],
                                                                  out_model['deletion_length_weights'])

    return out_model


def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before


class SequencingErrorModel:
    """
    This class represents the parameters needed to model quality scores.
    """
    def __init__(self, options):
        """
        This is a function version of the old SequencingErrors class. May convert to a class
        if it seems valuable. For now, let's try this.

        New input model:

        Dictionary of the probabilities and values:
            - quality_scores: list of possible values of quality scores.
            - quality_score_probabilities: 2D numpy list, with each row representing the
                                           probabilities of each quality bin at that location.
            - quality_offset: int representing the offset used to translate from quality score to ascii character,
                              usually this is 33, but some machines differ.b
            - read_length: int of the observed read length from the data.
            - average_error: float of the estimated average error of the entire dataset
            - is_uniform: a boolean that indicates if the quality scores have uniform error.
            - error_parameters:
                - error_params[0] - A nucleotide substitution matrix
                - error_params[1] - The chance that he error is an indel
                - error_params[2] - The probability distribution of the possible lengths of indels
                - error_params[3] - The possible lengths of indels
                - error_params[4] - The probability of an indel being an insertion
                - error_params[5] - Probability distribution for the 4 nucleotides for insertions

        """
        mssg = "Problem loading Sequencing error model data @error_model"
        error_model = pickle_load_model(gzip.open(options.error_model, 'rb'), mssg)
        # Convert qual scores to numpy array for faster recall.
        self.quality_scores = np.array(error_model['quality_scores'])
        # pre-compute the error rate for each quality score
        self.quality_score_error_rate = {x: 10. ** (-x / 10) for x in error_model['quality_scores']}
        self.average_error = error_model['average_error']
        self.read_length = error_model['read_length']
        self.uniform_quality_score = None
        self.rescale_qualities = options.rescale_qualities
        if error_model['is_uniform']:
            # Near as I can tell, this number is the average error translated into a phred score, but plus 0.5
            # for unknown reasons (legacy code). After running some tests, the effect is that sometimes the 0.5 bumps it
            # up, so I think this is basically a way to round it.
            avg_error_phred_score = int(-10. * np.log10(self.average_error) + 0.5)
            phred_score_bin = take_closest(self.quality_scores, avg_error_phred_score)
            self.uniform_quality_score = min([max(self.quality_scores), phred_score_bin])
        quality_score_probability_matrix = error_model['quality_score_probabilities']
        self.quality_score_probabilities = quality_score_probability_matrix.apply(
                lambda row: DiscreteDistribution(self.quality_scores, row), axis=1).to_numpy()
        matrix = pd.DataFrame(error_model['error_parameters'][0])
        self.substitution_model = matrix.apply(lambda row: DiscreteDistribution(ALLOWED_NUCL, row),
                                               axis=1).to_numpy()
        # Probability that if there is an error, it is an indel (v snp)
        self.indel_probability = error_model['error_parameters'][1]
        # Probability that an indel is an insertion
        self.insertion_probability = error_model['error_parameters'][4]
        # Possible lengths for indels, along with the weights of each
        self.indel_lengths = np.array(error_model['error_parameters'][3])
        self.indel_length_model = DiscreteDistribution(self.indel_lengths,
                                                       np.array(error_model['error_parameters'][2]))
        # Nucleotides and their probability values for a random choice
        self.nucleotide_insertion_model = error_model['error_parameters'][5]
        self.quality_offset = error_model['quality_offset']

        # Not sure what effect rescaling the quality scores will have on binned qualities.
        self.error_scale = 1.0
        if options.rescale_qualities:
            self.error_scale = options.avg_seq_error / error_model['average_error']

        # The point of this is supposed to be spread out the error model across a broader read length,
        # if the read length the user requests is different from what the sequencing error model calculated.
        self.quality_index_remap = range(options.read_len)
        if options.read_len != len(error_model['quality_score_probabilities']):
            self.quality_index_remap = np.array([max([0, len(error_model["quality_score_probabilities"]) * n //
                                                options.read_len]) for n in range(options.read_len)])

    def get_sequencing_errors(self, read_data, is_reverse_strand=False):
        """
        Inserts errors of type substitution, insertion, or deletion into read_data, and assigns a quality score
        based on the container model.
        :param read_data: sequence to insert errors into
        :param is_reverse_strand: whether to treat this as the reverse strand or not
        :return: modified sequence and associated quality scores
        """
        out_qualities = []
        introduced_errors = []
        sequencing_errors = []
        if self.uniform_quality_score:
            # Since this is uniform, forward and reverse are irrelevant
            out_qualities = [chr(self.uniform_quality_score + self.quality_offset)] * self.read_length
            # Let's only bother to do this if they even want errors. If error_scale was set to 0, then skip
            if self.error_scale != 0:
                for i in range(self.read_length):
                    if random.random() < self.error_scale * self.quality_score_error_rate[self.uniform_quality_score]:
                        sequencing_errors.append(i)
        else:
            temp_quality_array = [self.quality_score_probabilities[n].sample()
                                  for n in range(len(self.quality_score_probabilities))]

            # Now we remap the quality scores to the array based on the remap defined above.
            for i in range(self.read_length):
                out_qualities.append(int(temp_quality_array[self.quality_index_remap[i]]))
            # If this is the reverse strand, we just flip the quality array
            if is_reverse_strand:
                out_qualities = out_qualities[::-1]

            # Now we find where the errors are. This part works whether the strand is forward or reverse.
            # Let's only bother to do this if they even want errors. If error_scale was set to 0, then skip
            if self.error_scale != 0:
                for i in range(self.read_length):
                    if random.random() < self.error_scale * self.quality_score_error_rate[out_qualities[i]]:
                        sequencing_errors.append(i)

            # We'll see if this has any effect on bins. I doubt it. But since this method allows for bins of size 1,
            # then it might work in that situation.
            if self.rescale_qualities:
                # First we rescale the quality score literals then convert back to phred score (with the 0.5
                # causing borderline cases to take the next highest number).
                out_qualities = [max([0, int(-10. * np.log10(self.error_scale * self.quality_score_error_rate[n]) + 0.5)]) for n in out_qualities]
                # Now rebin the quality scores.
                out_qualities = [take_closest(self.quality_scores, n) for n in out_qualities]
                out_qualities = ''.join([chr(n + self.quality_offset) for n in out_qualities])
            else:
                out_qualities = ''.join([chr(n + self.quality_offset) for n in out_qualities])

        # The use case here would be someone running a simulation where they want no sequencing errors.
        # Since we skipped making errors if the error_scale was 0, sequencing_errors will be empty in this use case.
        if self.error_scale == 0:
            return out_qualities, sequencing_errors

        number_deletions_so_far = 0
        # don't allow indel errors to occur on subsequent positions
        previous_indel = -2
        # don't allow other sequencing errors to occur on bases removed by deletion errors
        del_blacklist = []

        for index in sequencing_errors[::-1]:
            # determine error type
            is_sub = True
            # This check checks that we are not trying to insert indels at the end of a read,
            # or overlapping with a different indel.
            if index != 0 and \
                    index != self.read_length - 1 - max(self.indel_lengths) and \
                    abs(index - previous_indel) > 1:
                if random.random() < self.indel_probability:
                    is_sub = False

            # Insert substitution error
            if is_sub:
                current_nucleotide = read_data[index]
                nuc_index = ALLOWED_NUCL.index(current_nucleotide)
                # take the zero index because this returns a list of length 1.
                new_nucleotide = self.substitution_model[nuc_index].sample()
                introduced_errors.append(('S', 1, index, current_nucleotide, new_nucleotide))

            # insert indel error:
            else:
                # Need to take the first element because this returns a list with 1 element.
                indel_len = self.indel_length_model.sample()

                # insertion error:
                if random.random() < self.insertion_probability:
                    current_nucleotide = read_data[index]
                    insert = random.choices(ALLOWED_NUCL, self.nucleotide_insertion_model, k=indel_len)
                    new_nucleotide = current_nucleotide + "".join(insert)
                    introduced_errors.append(('I', len(new_nucleotide) - 1, index, current_nucleotide, new_nucleotide))

                elif index < self.read_length - 2 - number_deletions_so_far:
                    current_nucleotide = read_data[index: index + indel_len + 1]
                    new_nucleotide = read_data[index]
                    number_deletions_so_far += len(current_nucleotide) - 1
                    introduced_errors.append(('D', len(current_nucleotide) - 1, index,
                                              current_nucleotide, new_nucleotide))
                    del_blacklist.extend(list(range(index + 1, index + indel_len + 1)))

                previous_indel = index

        # Remove blacklisted errors
        for i in range(len(introduced_errors) - 1, -1, -1):
            if introduced_errors[i][2] in del_blacklist:
                del introduced_errors[i]

        return out_qualities, introduced_errors


class Models:
    """
    Sets up models for use in the rest of the program.
    """
    def __init__(self, options):

        # Load mutation models
        self.mutation_model = parse_input_mutation_model(options.mutation_model)

        self.cancer_model = None
        if options.cancer:
            self.cancer_model = parse_input_mutation_model(options.cancer_model, True)

        log_mssg("Mutation models loaded", 'debug')

        # We need sequencing errors to get the quality score attributes, even for the vcf
        self.sequencing_error_model = SequencingErrorModel(options)

        log_mssg('Sequencing error model loaded', 'debug')

        mssg = "f'ERROR: problem reading @gc_model. Please check file path and try again. " \
               "This file should be the output of compute_gc.py'"
        self.gc_model = pickle_load_model(gzip.open(options.gc_model, 'r'), mssg)

        log_mssg('GC Bias model loaded', 'debug')

        self.fraglen_model = None

        if options.paired_ended:
            if options.fragment_model:
                mssg = 'Problem loading the empirical fragment length model @fragment_model. Please check file and try' \
                       'again.'
                log_mssg("Using empirical fragment length distribution", 'info')
                potential_values, potential_prob = pickle_load_model(open(options.fragment_model, 'rb'), mssg)

                fraglen_values = []
                fraglen_probability = []
                for i in range(len(potential_values)):
                    if potential_values[1] > options.read_len:
                        fraglen_values.append(potential_values[i])
                        fraglen_probability.append(potential_prob[i])

                self.fraglen_model = DiscreteDistribution(fraglen_values, fraglen_probability)
                options.set_value('fragment_mean', fraglen_values[mean_ind_of_weighted_list(fraglen_probability)])

            # Using artificial fragment length distribution, if the parameters were specified
            # fragment length distribution: normal distribution that goes out to +- 6 standard deviations
            else:
                log_mssg(f'Using artificial fragment length distribution.', 'info')
                if options.fragment_st_dev == 0:
                    self.fraglen_model = DiscreteDistribution([options.fragment_mean], [1],
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
                    self.fraglen_model = DiscreteDistribution(fraglen_values, fraglen_probability)
            log_mssg(f'Loaded paired-end models', 'debug')
