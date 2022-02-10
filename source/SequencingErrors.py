import numpy as np
import random
import pickle
from source.probability import DiscreteDistribution
from source.constants_and_models import ALLOWED_NUCL, NUC_IND
from source.error_handling import premature_exit, print_and_log


class SequencingErrors:
    """
    Container for read data: computes quality scores and positions to insert errors
    """

    def __init__(self, read_len, error_model, rescaled_error, rescale_qual=False, debug=False):

        self.read_len = read_len
        self.rescale_qual = rescale_qual

        self.uniform = False

        try:
            error_dat = pickle.load(open(error_model, 'rb'), encoding="bytes")
        except IOError:
            print("\nProblem opening the sequencing error model.\n")
            sys.exit(1)

        self.uniform = False

        # uniform-error SE reads (e.g., PacBio)
        if len(error_dat) == 4:
            self.uniform = True
            [q_scores, off_q, avg_error, error_params] = error_dat
            self.uniform_q_score = min([max(q_scores), int(-10. * np.log10(avg_error) + 0.5)])
            print_and_log(f'Reading in uniform sequencing error model... (q={self.uniform_q_score}+{off_q}, '
                          f'p(err)={(100. * avg_error):0.2f}%)', 'info')

        # only 1 q-score model present, use same model for both strands
        elif len(error_dat) == 6:
            [init_q1, prob_q1, q_scores, off_q, avg_error, error_params] = error_dat
            self.pe_models = False

        # found a q-score model for both forward and reverse strands
        elif len(error_dat) == 8:
            [init_q1, prob_q1, init_q2, prob_q2, q_scores, off_q, avg_error, error_params] = error_dat
            self.pe_models = True
            if len(init_q1) != len(init_q2) or len(prob_q1) != len(prob_q2):
                print_and_log(f'R1 and R2 quality score models are of different length.', 'error')
                premature_exit(1)

        # This serves as a sanity check for the input model
        else:
            print_and_log('Something wrong with error model.', 'error')
            if debug:
                print_and_log(f"error model had a length of {len(error_model)}", 'debug')
            premature_exit(1)

        self.q_err_rate = [0.] * (max(q_scores) + 1)
        for q in q_scores:
            self.q_err_rate[q] = 10. ** (-q / 10.)
        self.off_q = off_q
        self.err_p = error_params

        # Selects a new nucleotide based on the error model
        self.err_sse = [DiscreteDistribution(n, ALLOWED_NUCL) for n in self.err_p[0]]

        # allows for selection of indel length based on the parameters of the model
        self.err_sie = DiscreteDistribution(self.err_p[2], self.err_p[3])

        # allows for indel insertion based on the length above and the probability from the model
        self.err_sin = DiscreteDistribution(self.err_p[5], ALLOWED_NUCL)

        # adjust sequencing error frequency to match desired rate
        if rescaled_error is None:
            self.error_scale = 1.0
        else:
            self.error_scale = rescaled_error / avg_error
            if not self.rescale_qual:
                print_and_log(f'Quality scores no longer exactly representative of error probability. '
                              f'Error model scaled by {self.error_scale:.3f} to match desired rate...', 'warning')
            if self.uniform:
                if rescaled_error <= 0.:
                    self.uniform_q_score = max(q_scores)
                else:
                    self.uniform_q_score = min([max(q_scores), int(-10. * np.log10(rescaled_error) + 0.5)])
                print_and_log(f'Uniform quality score scaled to match specified error rate '
                              f'(q={self.uniform_qscore}+{self.off_q}, p(err)={(100. * rescaled_error):0.2f}%)', 'info')

        if not self.uniform:
            # adjust length to match desired read length
            if self.read_len == len(init_q1):
                self.q_ind_remap = range(self.read_len)
            else:
                print_and_log(f'Read length of error model ({len(init_q1)}) '
                              f'does not match -R value ({self.read_len}), rescaling model...', 'warning')
                self.q_ind_remap = [max([1, len(init_q1) * n // read_len]) for n in range(read_len)]

            # initialize probability distributions
            self.init_dist_by_pos_1 = [DiscreteDistribution(init_q1[i], q_scores) for i in range(len(init_q1))]
            self.prob_dist_by_pos_by_prev_q1 = [None]
            for i in range(1, len(init_q1)):
                self.prob_dist_by_pos_by_prev_q1.append([])
                for j in range(len(init_q1[0])):
                    # if we don't have sufficient data for a transition, use the previous quality score
                    if np.sum(prob_q1[i][j]) <= 0.:
                        self.prob_dist_by_pos_by_prev_q1[-1].append(
                            DiscreteDistribution([1], [q_scores[j]], degenerate_val=q_scores[j]))
                    else:
                        self.prob_dist_by_pos_by_prev_q1[-1].append(DiscreteDistribution(prob_q1[i][j], q_scores))

            # If paired-end, initialize probability distributions for the other strand
            if self.pe_models:
                self.init_dist_by_pos_2 = [DiscreteDistribution(init_q2[i], q_scores) for i in range(len(init_q2))]
                self.prob_dist_by_pos_by_prev_q2 = [None]
                for i in range(1, len(init_q2)):
                    self.prob_dist_by_pos_by_prev_q2.append([])
                    for j in range(len(init_q2[0])):
                        # if we don't have sufficient data for a transition, use the previous qscore
                        if np.sum(prob_q2[i][j]) <= 0.:
                            self.prob_dist_by_pos_by_prev_q2[-1].append(
                                DiscreteDistribution([1], [q_scores[j]], degenerate_val=q_scores[j]))
                        else:
                            self.prob_dist_by_pos_by_prev_q2[-1].append(DiscreteDistribution(prob_q2[i][j], q_scores))

    def get_sequencing_errors(self, read_data, is_reverse_strand=False):
        """
        Inserts errors of type substitution, insertion, or deletion into read_data, and assigns a quality score
        based on the container model.
        :param read_data: sequence to insert errors into
        :param is_reverse_strand: whether to treat this as the reverse strand or not
        :return: modified sequence and associate quality scores
        """

        # TODO this is one of the slowest methods in the code. Need to investigate how to speed this up.
        q_out = [0] * self.read_len
        s_err = []

        if self.uniform:
            my_q = [self.uniform_q_score + self.off_q] * self.read_len
            q_out = ''.join([chr(n) for n in my_q])
            for i in range(self.read_len):
                if random.random() < self.error_scale * self.q_err_rate[self.uniform_q_score]:
                    s_err.append(i)
        else:
            if self.pe_models and is_reverse_strand:
                my_q = self.init_dist_by_pos_2[0].sample()
            else:
                my_q = self.init_dist_by_pos_1[0].sample()
            q_out[0] = my_q

            # Every time this is hit, we loop the entire read length twice. I feel like these two loops
            # Could be combined into one fairly easily. The culprit seems to bee too many hits to the sample() method.
            for i in range(1, self.read_len):
                if self.pe_models and is_reverse_strand:
                    my_q = self.prob_dist_by_pos_by_prev_q2[self.q_ind_remap[i]][my_q].sample()
                else:
                    my_q = self.prob_dist_by_pos_by_prev_q1[self.q_ind_remap[i]][my_q].sample()
                q_out[i] = my_q

            if is_reverse_strand:
                q_out = q_out[::-1]

            for i in range(self.read_len):
                if random.random() < self.error_scale * self.q_err_rate[q_out[i]]:
                    s_err.append(i)

            if self.rescale_qual:  # do we want to rescale qual scores to match rescaled error?
                q_out = [max([0, int(-10. * np.log10(self.error_scale * self.q_err_rate[n]) + 0.5)]) for n in q_out]
                q_out = [min([int(self.q_err_rate[-1]), n]) for n in q_out]
                q_out = ''.join([chr(n + self.off_q) for n in q_out])
            else:
                q_out = ''.join([chr(n + self.off_q) for n in q_out])

        if self.error_scale == 0.0:
            return q_out, []

        s_out = []
        n_del_so_far = 0
        # don't allow indel errors to occur on subsequent positions
        prev_indel = -2
        # don't allow other sequencing errors to occur on bases removed by deletion errors
        del_blacklist = []

        # Need to check into this loop, to make sure it isn't slowing us down.
        # The culprit seems to bee too many hits to the sample() method. This has a few of those calls.
        for ind in s_err[::-1]:  # for each error that we're going to insert...

            # determine error type
            is_sub = True
            if ind != 0 and ind != self.read_len - 1 - max(self.err_p[3]) and abs(ind - prev_indel) > 1:
                if random.random() < self.err_p[1]:
                    is_sub = False

            # insert substitution error
            if is_sub:
                my_nucl = str(read_data[ind])
                new_nucl = self.err_sse[NUC_IND[my_nucl]].sample()
                s_out.append(('S', 1, ind, my_nucl, new_nucl))

            # insert indel error
            else:
                indel_len = self.err_sie.sample()

                # insertion error
                if random.random() < self.err_p[4]:
                    my_nucl = str(read_data[ind])
                    new_nucl = my_nucl + ''.join([self.err_sin.sample() for n in range(indel_len)])
                    s_out.append(('I', len(new_nucl) - 1, ind, my_nucl, new_nucl))

                # deletion error (prevent too many of them from stacking up)
                elif ind < self.read_len - 2 - n_del_so_far:
                    my_nucl = str(read_data[ind:ind + indel_len + 1])
                    new_nucl = str(read_data[ind])
                    n_del_so_far += len(my_nucl) - 1
                    s_out.append(('D', len(my_nucl) - 1, ind, my_nucl, new_nucl))
                    for i in range(ind + 1, ind + indel_len + 1):
                        del_blacklist.append(i)
                prev_indel = ind

        # remove blacklisted errors
        for i in range(len(s_out) - 1, -1, -1):
            if s_out[i][2] in del_blacklist:
                del s_out[i]

        return q_out, s_out
