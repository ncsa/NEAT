import random
import copy
import bisect
import sys

import numpy as np
from Bio.Seq import Seq
from Bio.Seq import MutableSeq

from source.neat_cigar import CigarString
from source.probability import DiscreteDistribution, poisson_list
from source.constants_and_models import ALLOWED_NUCL, DEFAULT_MODEL_1, ALL_IND, COV_FRAGLEN_PERCENTILE, TRI_IND
from source.constants_and_models import LARGE_NUMBER, MAX_MUTFRAC, IGNORE_TRINUC, MAX_ATTEMPTS, NUC_IND


# TODO This whole file is in desperate need of refactoring
# TODO 1st step will be to make all cigar strings optional, since that is only needed for bams


class SequenceContainer:
    """
    Container for reference sequences, applies mutations
    """

    def __init__(self, x_offset, sequence, ploidy, window_overlap, read_len, mut_models=None,
                 mut_rate=None, only_vcf=False):

        # initialize basic variables
        self.only_vcf = only_vcf
        self.x = x_offset
        self.ploidy = ploidy
        self.read_len = read_len
        self.sequences = [sequence] * self.ploidy
        self.seq_len = len(sequence)
        self.indel_list = [[]] * self.ploidy
        self.snp_list = [[]] * self.ploidy
        self.fm_pos = [[]] * self.ploidy
        self.fm_span = [[]] * self.ploidy
        self.all_cigar = [[]] * self.ploidy

        # Blacklist explanation:
        # black_list[ploid][pos] = 0		safe to insert variant here
        # black_list[ploid][pos] = 1		indel inserted here
        # black_list[ploid][pos] = 2		snp inserted here
        # black_list[ploid][pos] = 3		invalid position for various processing reasons
        self.black_list = [np.zeros(self.seq_len, dtype='<i4')] * self.ploidy

        # disallow mutations to occur on window overlap points
        self.win_buffer = window_overlap
        for p in range(self.ploidy):
            self.black_list[p][-self.win_buffer] = 3
            self.black_list[p][-self.win_buffer - 1] = 3

        # initialize mutation models
        if not mut_models:
            default_model = [copy.deepcopy(DEFAULT_MODEL_1)] * self.ploidy
            self.model_data = default_model[:self.ploidy]
        else:
            if len(mut_models) != self.ploidy:
                print('\nError: Number of mutation models received is not equal to specified ploidy\n')
                sys.exit(1)
            self.model_data = copy.deepcopy(mut_models)

        # do we need to rescale mutation frequencies?
        mut_rate_sum = sum([n[0] for n in self.model_data])
        self.mut_rescale = mut_rate
        if self.mut_rescale is None:
            self.mut_scalar = 1.0
        else:
            self.mut_scalar = float(self.mut_rescale) / (mut_rate_sum / float(len(self.model_data)))

        # how are mutations spread to each ploid, based on their specified mut rates?
        self.ploid_mut_frac = [float(n[0]) / mut_rate_sum for n in self.model_data]
        self.ploid_mut_prior = DiscreteDistribution(self.ploid_mut_frac, range(self.ploidy))

        # init mutation models
        #
        # self.models[ploid][0] = average mutation rate
        # self.models[ploid][1] = p(mut is homozygous | mutation occurs)
        # self.models[ploid][2] = p(mut is indel | mut occurs)
        # self.models[ploid][3] = p(insertion | indel occurs)
        # self.models[ploid][4] = distribution of insertion lengths
        # self.models[ploid][5] = distribution of deletion lengths
        # self.models[ploid][6] = distribution of trinucleotide SNP transitions
        # self.models[ploid][7] = p(trinuc mutates)
        self.models = []
        for n in self.model_data:
            self.models.append([self.mut_scalar * n[0], n[1], n[2], n[3], DiscreteDistribution(n[5], n[4]),
                                DiscreteDistribution(n[7], n[6]), []])
            for m in n[8]:
                # noinspection PyTypeChecker
                self.models[-1][6].append(
                    [DiscreteDistribution(m[0], ALLOWED_NUCL), DiscreteDistribution(m[1], ALLOWED_NUCL),
                     DiscreteDistribution(m[2], ALLOWED_NUCL), DiscreteDistribution(m[3], ALLOWED_NUCL)])
            self.models[-1].append([m for m in n[9]])

        # initialize poisson attributes
        self.indel_poisson, self.snp_poisson = self.init_poisson()

        # sample the number of variants that will be inserted into each ploid
        self.indels_to_add = [n.sample() for n in self.indel_poisson]
        self.snps_to_add = [n.sample() for n in self.snp_poisson]

        # initialize trinuc snp bias
        # compute mutation positional bias given trinucleotide strings of the sequence (ONLY AFFECTS SNPs)
        #
        # note: since indels are added before snps, it's possible these positional biases aren't correctly utilized
        #       at positions affected by indels. At the moment I'm going to consider this negligible.
        trinuc_snp_bias = [[0.] * self.seq_len] * self.ploidy
        self.trinuc_bias = [None] * self.ploidy
        for p in range(self.ploidy):
            for i in range(self.win_buffer + 1, self.seq_len - 1):
                trinuc_snp_bias[p][i] = self.models[p][7][ALL_IND[str(self.sequences[p][i - 1:i + 2])]]
            self.trinuc_bias[p] = DiscreteDistribution(trinuc_snp_bias[p][self.win_buffer + 1:self.seq_len - 1],
                                                       range(self.win_buffer + 1, self.seq_len - 1))

        # initialize coverage attributes
        self.window_size = None
        self.coverage_distribution = None
        self.fraglen_ind_map = None

    def update_basic_vars(self, x_offset, sequence, ploidy, window_overlap, read_len):
        self.x = x_offset
        self.ploidy = ploidy
        self.read_len = read_len
        self.sequences = [sequence] * self.ploidy
        self.seq_len = len(sequence)
        self.indel_list = [[]] * self.ploidy
        self.snp_list = [[]] * self.ploidy
        self.fm_pos = [[]] * self.ploidy
        self.fm_span = [[]] * self.ploidy
        self.black_list = [np.zeros(self.seq_len, dtype='<i4')] * self.ploidy
        self.all_cigar = [[]] * self.ploidy

        # disallow mutations to occur on window overlap points
        self.win_buffer = window_overlap
        for p in range(self.ploidy):
            self.black_list[p][-self.win_buffer] = 3
            self.black_list[p][-self.win_buffer - 1] = 3

    def update_mut_models(self, mut_models, mut_rate):
        if not mut_models:
            default_model = [copy.deepcopy(DEFAULT_MODEL_1)] * self.ploidy
            self.model_data = default_model[:self.ploidy]
        else:
            if len(mut_models) != self.ploidy:
                print('\nError: Number of mutation models received is not equal to specified ploidy\n')
                sys.exit(1)
            self.model_data = copy.deepcopy(mut_models)

        # do we need to rescale mutation frequencies?
        mut_rate_sum = sum([n[0] for n in self.model_data])
        self.mut_rescale = mut_rate
        if self.mut_rescale is None:
            self.mut_scalar = 1.0
        else:
            self.mut_scalar = float(self.mut_rescale) / (mut_rate_sum / float(len(self.model_data)))

        # how are mutations spread to each ploid, based on their specified mut rates?
        self.ploid_mut_frac = [float(n[0]) / mut_rate_sum for n in self.model_data]
        self.ploid_mut_prior = DiscreteDistribution(self.ploid_mut_frac, range(self.ploidy))
        self.models = []
        for n in self.model_data:
            self.models.append([self.mut_scalar * n[0], n[1], n[2], n[3], DiscreteDistribution(n[5], n[4]),
                                DiscreteDistribution(n[7], n[6]), []])
            for m in n[8]:
                # noinspection PyTypeChecker
                self.models[-1][6].append([DiscreteDistribution(m[0], ALLOWED_NUCL),
                                           DiscreteDistribution(m[1], ALLOWED_NUCL),
                                           DiscreteDistribution(m[2], ALLOWED_NUCL),
                                           DiscreteDistribution(m[3], ALLOWED_NUCL)])
            self.models[-1].append([m for m in n[9]])

    def update_trinuc_bias(self):
        trinuc_snp_bias = [[0.] * self.seq_len] * self.ploidy
        self.trinuc_bias = [None] * self.ploidy
        for p in range(self.ploidy):
            for i in range(self.win_buffer + 1, self.seq_len - 1):
                trinuc_snp_bias[p][i] = self.models[p][7][ALL_IND[str(self.sequences[p][i - 1:i + 2])]]
            self.trinuc_bias[p] = DiscreteDistribution(trinuc_snp_bias[p][self.win_buffer + 1:self.seq_len - 1],
                                                       range(self.win_buffer + 1, self.seq_len - 1))

    def init_coverage_new(self, coverage_data, frag_dist=None):
        """
        Initializes coverage for the sequence container. Only makes changes if we are not in vcf-only mode.

        :param coverage_data: A tuple containing the window size, gc scalars and target coverage values.
        :param frag_dist: A probability distribution of the fragment size.
        :return: Mean coverage value
        """
        # TODO This new form of init_coverage is attempting to divorce the process from the cigar string creation,
        #  but will require some further investigation.

        # TODO this section is also quite slow
        # If we're only creating a vcf, skip some expensive initialization related to coverage depth
        if not self.only_vcf:
            (self.window_size, gc_scalars, target_cov_vals) = coverage_data
            gc_cov_vals = [[]] * self.ploidy
            tr_cov_vals = [[]] * self.ploidy
            avg_out = []
            self.coverage_distribution = []
            for i in range(self.ploidy):
                # Combined a couple of lines for this. I'm trying to divorce these calculations from the cigar creation
                if len(self.sequences[i]) - self.read_len > 0:
                    max_coord = len(self.sequences[i]) - self.read_len
                else:
                    max_coord = len(self.sequences[i])

                # compute gc-bias
                j = 0
                while j + self.window_size < len(self.sequences[i]):
                    gc_c = self.sequences[i][j:j + self.window_size].count('G') + \
                           self.sequences[i][j:j + self.window_size].count('C')
                    gc_cov_vals[i].extend([gc_scalars[gc_c]] * self.window_size)
                    j += self.window_size
                gc_c = self.sequences[i][-self.window_size:].count('G') + \
                       self.sequences[i][-self.window_size:].count('C')
                gc_cov_vals[i].extend([gc_scalars[gc_c]] * (len(self.sequences[i]) - len(gc_cov_vals[i])))

                # Targeted values
                tr_cov_vals[i].append(target_cov_vals[0])

                prev_val = self.fm_pos[i][0]
                for j in range(1, max_coord):
                    if self.fm_pos[i][j] is None:
                        tr_cov_vals[i].append(target_cov_vals[prev_val])
                    elif self.fm_span[i][j] - self.fm_pos[i][j] <= 1:
                        tr_cov_vals[i].append(target_cov_vals[prev_val])
                    else:
                        tr_cov_vals[i].append(sum(target_cov_vals[self.fm_pos[i][j]:self.fm_span[i][j]]) / float(
                            self.fm_span[i][j] - self.fm_pos[i][j]))
                        prev_val = self.fm_pos[i][j]
                    # TODO add debug attribute to this class to activate statements like this
                    # Debug statement
                    # print(f'({i, j}), {self.all_cigar[i][j]}, {self.fm_pos[i][j]}, {self.fm_span[i][j]}')

                # shift by half of read length
                if len(tr_cov_vals[i]) > int(self.read_len / 2.):
                    tr_cov_vals[i] = [0.0] * int(self.read_len // 2) + tr_cov_vals[i][:-int(self.read_len / 2.)]
                # fill in missing indices
                tr_cov_vals[i].extend([0.0] * (len(self.sequences[i]) - len(tr_cov_vals[i])))

                #
                coverage_vector = np.cumsum([tr_cov_vals[i][nnn] *
                                             gc_cov_vals[i][nnn] for nnn in range(len(tr_cov_vals[i]))])
                coverage_vals = []
                # TODO if max_coord is <=0, this is a problem
                for j in range(0, max_coord):
                    coverage_vals.append(coverage_vector[j + self.read_len] - coverage_vector[j])
                # Below is Zach's attempt to fix this. The commented out line is the original
                # avg_out.append(np.mean(coverage_vals) / float(self.read_len))
                avg_out.append(np.mean(coverage_vals) / float(min([self.read_len, max_coord])))
                # Debug statement
                # print(f'{avg_out}, {np.mean(avg_out)}')

                if frag_dist is None:
                    # Debug statement
                    # print(f'++++, {max_coord}, {len(self.sequences[i])}, '
                    #       f'{len(self.all_cigar[i])}, {len(coverage_vals)}')
                    self.coverage_distribution.append(DiscreteDistribution(coverage_vals, range(len(coverage_vals))))

                # fragment length nightmare
                else:
                    current_thresh = 0.
                    index_list = [0]
                    for j in range(len(frag_dist.cum_prob)):
                        if frag_dist.cum_prob[j] >= current_thresh + COV_FRAGLEN_PERCENTILE / 100.0:
                            current_thresh = frag_dist.cum_prob[j]
                            index_list.append(j)
                    flq = [frag_dist.values[nnn] for nnn in index_list]
                    if frag_dist.values[-1] not in flq:
                        flq.append(frag_dist.values[-1])
                    flq.append(LARGE_NUMBER)

                    self.fraglen_ind_map = {}
                    for j in frag_dist.values:
                        b_ind = bisect.bisect(flq, j)
                        if abs(flq[b_ind - 1] - j) <= abs(flq[b_ind] - j):
                            self.fraglen_ind_map[j] = flq[b_ind - 1]
                        else:
                            self.fraglen_ind_map[j] = flq[b_ind]

                    self.coverage_distribution.append({})
                    for flv in sorted(list(set(self.fraglen_ind_map.values()))):
                        buffer_val = self.read_len
                        for j in frag_dist.values:
                            if self.fraglen_ind_map[j] == flv and j > buffer_val:
                                buffer_val = j
                        max_coord = min([len(self.sequences[i]) - buffer_val - 1,
                                         len(self.all_cigar[i]) - buffer_val + self.read_len - 2])
                        # print 'BEFORE:', len(self.sequences[i])-buffer_val
                        # print 'AFTER: ', len(self.all_cigar[i])-buffer_val+self.read_len-2
                        # print 'AFTER2:', max_coord
                        coverage_vals = []
                        for j in range(0, max_coord):
                            coverage_vals.append(
                                coverage_vector[j + self.read_len] - coverage_vector[j] + coverage_vector[j + flv] -
                                coverage_vector[
                                    j + flv - self.read_len])

                        # EXPERIMENTAL
                        # quantized_cov_vals = quantize_list(coverage_vals)
                        # self.coverage_distribution[i][flv] = \
                        #     DiscreteDistribution([n[2] for n in quantized_cov_vals],
                        #                          [(n[0], n[1]) for n in quantized_cov_vals])

                        # TESTING
                        # import matplotlib.pyplot as mpl
                        # print len(coverage_vals),'-->',len(quantized_cov_vals)
                        # mpl.figure(0)
                        # mpl.plot(range(len(coverage_vals)), coverage_vals)
                        # for qcv in quantized_cov_vals:
                        # mpl.plot([qcv[0], qcv[1]+1], [qcv[2],qcv[2]], 'r')
                        # mpl.show()
                        # sys.exit(1)

                        self.coverage_distribution[i][flv] = DiscreteDistribution(coverage_vals,
                                                                                  range(len(coverage_vals)))

            return np.mean(avg_out)

    def init_coverage(self, coverage_data, frag_dist=None):
        """
        Initializes coverage for the sequence container. Only makes changes if we are not in vcf-only mode.
        :param coverage_data: A tuple containing the window size, gc scalars and target coverage values.
        :param frag_dist: A probability distribution of the fragment size.
        :return: Mean coverage value
        """

        # TODO this section is also quite slow and will need further investigation
        # If we're only creating a vcf, skip some expensive initialization related to coverage depth
        if not self.only_vcf:
            (self.window_size, gc_scalars, target_cov_vals) = coverage_data
            gc_cov_vals = [[] for _ in self.sequences]
            tr_cov_vals = [[] for _ in self.sequences]
            avg_out = []
            self.coverage_distribution = []
            for i in range(len(self.sequences)):
                # Zach implemented a change here but I can't remember if I changed it back for some reason.
                # If second line below doesn't work, reactivate the first line.
                # max_coord = min([len(self.sequences[i]) - self.read_len, len(self.all_cigar[i]) - self.read_len])
                max_coord = min([len(self.sequences[i]) - self.read_len, len(self.all_cigar[i]) - 1])

                # Trying to fix a problem wherein the above line gives a negative answer
                if max_coord <= 0:
                    max_coord = min([len(self.sequences[i]), len(self.all_cigar[i])])

                # compute gc-bias
                j = 0
                while j + self.window_size < len(self.sequences[i]):
                    gc_c = self.sequences[i][j:j + self.window_size].count('G') + \
                           self.sequences[i][j:j + self.window_size].count('C')
                    gc_cov_vals[i].extend([gc_scalars[gc_c]] * self.window_size)
                    j += self.window_size
                gc_c = self.sequences[i][-self.window_size:].count('G') + \
                       self.sequences[i][-self.window_size:].count('C')
                gc_cov_vals[i].extend([gc_scalars[gc_c]] * (len(self.sequences[i]) - len(gc_cov_vals[i])))

                # Targeted values
                tr_cov_vals[i].append(target_cov_vals[0])
                prev_val = self.fm_pos[i][0]
                for j in range(1, max_coord):
                    if self.fm_pos[i][j] is None:
                        tr_cov_vals[i].append(target_cov_vals[prev_val])
                    elif self.fm_span[i][j] - self.fm_pos[i][j] <= 1:
                        tr_cov_vals[i].append(target_cov_vals[prev_val])
                    else:
                        tr_cov_vals[i].append(sum(target_cov_vals[self.fm_pos[i][j]:self.fm_span[i][j]]) / float(
                            self.fm_span[i][j] - self.fm_pos[i][j]))
                        prev_val = self.fm_pos[i][j]
                    # Debug statement
                    # print(f'({i, j}), {self.all_cigar[i][j]}, {self.fm_pos[i][j]}, {self.fm_span[i][j]}')

                # shift by half of read length
                if len(tr_cov_vals[i]) > int(self.read_len / 2.):
                    tr_cov_vals[i] = [0.0] * int(self.read_len // 2) + tr_cov_vals[i][:-int(self.read_len / 2.)]
                # fill in missing indices
                tr_cov_vals[i].extend([0.0] * (len(self.sequences[i]) - len(tr_cov_vals[i])))

                #
                coverage_vector = np.cumsum([tr_cov_vals[i][nnn] *
                                             gc_cov_vals[i][nnn] for nnn in range(len(tr_cov_vals[i]))])
                coverage_vals = []
                # TODO if max_coord is <=0, this is a problem
                for j in range(0, max_coord):
                    coverage_vals.append(coverage_vector[j + self.read_len] - coverage_vector[j])
                # Below is Zach's attempt to fix this. The commented out line is the original
                # avg_out.append(np.mean(coverage_vals) / float(self.read_len))
                avg_out.append(np.mean(coverage_vals) / float(min([self.read_len, max_coord])))
                # Debug statement
                # print(f'{avg_out}, {np.mean(avg_out)}')

                if frag_dist is None:
                    # Debug statement
                    # print(f'++++, {max_coord}, {len(self.sequences[i])}, '
                    #       f'{len(self.all_cigar[i])}, {len(coverage_vals)}')
                    self.coverage_distribution.append(DiscreteDistribution(coverage_vals, range(len(coverage_vals))))

                # fragment length nightmare
                else:
                    current_thresh = 0.
                    index_list = [0]
                    for j in range(len(frag_dist.cum_prob)):
                        if frag_dist.cum_prob[j] >= current_thresh + COV_FRAGLEN_PERCENTILE / 100.0:
                            current_thresh = frag_dist.cum_prob[j]
                            index_list.append(j)
                    flq = [frag_dist.values[nnn] for nnn in index_list]
                    if frag_dist.values[-1] not in flq:
                        flq.append(frag_dist.values[-1])
                    flq.append(LARGE_NUMBER)

                    self.fraglen_ind_map = {}
                    for j in frag_dist.values:
                        b_ind = bisect.bisect(flq, j)
                        if abs(flq[b_ind - 1] - j) <= abs(flq[b_ind] - j):
                            self.fraglen_ind_map[j] = flq[b_ind - 1]
                        else:
                            self.fraglen_ind_map[j] = flq[b_ind]

                    self.coverage_distribution.append({})
                    for flv in sorted(list(set(self.fraglen_ind_map.values()))):
                        buffer_val = self.read_len
                        for j in frag_dist.values:
                            if self.fraglen_ind_map[j] == flv and j > buffer_val:
                                buffer_val = j
                        max_coord = min([len(self.sequences[i]) - buffer_val - 1,
                                         len(self.all_cigar[i]) - buffer_val + self.read_len - 2])
                        # print 'BEFORE:', len(self.sequences[i])-buffer_val
                        # print 'AFTER: ', len(self.all_cigar[i])-buffer_val+self.read_len-2
                        # print 'AFTER2:', max_coord
                        coverage_vals = []
                        for j in range(0, max_coord):
                            coverage_vals.append(
                                coverage_vector[j + self.read_len] - coverage_vector[j] + coverage_vector[j + flv] -
                                coverage_vector[
                                    j + flv - self.read_len])

                        # EXPERIMENTAL
                        # quantized_cov_vals = quantize_list(coverage_vals)
                        # self.coverage_distribution[i][flv] = \
                        #     DiscreteDistribution([n[2] for n in quantized_cov_vals],
                        #                          [(n[0], n[1]) for n in quantized_cov_vals])

                        # TESTING
                        # import matplotlib.pyplot as mpl
                        # print len(coverage_vals),'-->',len(quantized_cov_vals)
                        # mpl.figure(0)
                        # mpl.plot(range(len(coverage_vals)), coverage_vals)
                        # for qcv in quantized_cov_vals:
                        # mpl.plot([qcv[0], qcv[1]+1], [qcv[2],qcv[2]], 'r')
                        # mpl.show()
                        # sys.exit(1)

                        self.coverage_distribution[i][flv] = DiscreteDistribution(coverage_vals,
                                                                                  range(len(coverage_vals)))

            return np.mean(avg_out)

    def init_poisson(self):
        ind_l_list = [self.seq_len * self.models[i][0] * self.models[i][2] * self.ploid_mut_frac[i] for i in
                      range(len(self.models))]
        snp_l_list = [self.seq_len * self.models[i][0] * (1. - self.models[i][2]) * self.ploid_mut_frac[i] for i in
                      range(len(self.models))]
        k_range = range(int(self.seq_len * MAX_MUTFRAC))
        # return (indel_poisson, snp_poisson)
        # TODO These next two lines are really slow. Maybe there's a better way
        return [poisson_list(k_range, ind_l_list[n]) for n in range(len(self.models))], \
               [poisson_list(k_range, snp_l_list[n]) for n in range(len(self.models))]

    def update(self, x_offset, sequence, ploidy, window_overlap, read_len, mut_models=None, mut_rate=None):
        # if mutation model is changed, we have to reinitialize it...
        if ploidy != self.ploidy or mut_rate != self.mut_rescale or mut_models is not None:
            self.ploidy = ploidy
            self.mut_rescale = mut_rate
            self.update_mut_models(mut_models, mut_rate)
        # if sequence length is different than previous window, we have to redo snp/indel poissons
        if len(sequence) != self.seq_len:
            self.seq_len = len(sequence)
            self.indel_poisson, self.snp_poisson = self.init_poisson()
        # basic vars
        self.update_basic_vars(x_offset, sequence, ploidy, window_overlap, read_len)
        self.indels_to_add = [n.sample() for n in self.indel_poisson]
        self.snps_to_add = [n.sample() for n in self.snp_poisson]
        # initialize trinuc snp bias
        if not IGNORE_TRINUC:
            self.update_trinuc_bias()

    def insert_mutations(self, input_list):
        for input_variable in input_list:
            which_ploid = []
            wps = input_variable[4][0]

            # if no genotype given, assume heterozygous and choose a single ploid based on their mut rates
            if wps is None:
                which_ploid.append(self.ploid_mut_prior.sample())
                which_alt = [0]
            else:
                if '/' in wps or '|' in wps:
                    if '/' in wps:
                        splt = wps.split('/')
                    else:
                        splt = wps.split('|')
                    which_ploid = []
                    for i in range(len(splt)):
                        if splt[i] == '1':
                            which_ploid.append(i)
                    # assume we're just using first alt for inserted variants?
                    which_alt = [0 for _ in which_ploid]
                # otherwise assume monoploidy
                else:
                    which_ploid = [0]
                    which_alt = [0]

            # ignore invalid ploids
            for i in range(len(which_ploid) - 1, -1, -1):
                if which_ploid[i] >= self.ploidy:
                    del which_ploid[i]

            for i in range(len(which_ploid)):
                p = which_ploid[i]
                my_alt = input_variable[2][which_alt[i]]
                my_var = (input_variable[0] - self.x, input_variable[1], my_alt)
                # This is a potential fix implemented by Zach in a previous commit. He left the next line in.
                # in_len = max([len(input_variable[1]), len(my_alt)])
                in_len = len(input_variable[1])

                if my_var[0] < 0 or my_var[0] >= len(self.black_list[p]):
                    print('\nError: Attempting to insert variant out of window bounds:')
                    print(my_var, '--> blackList[0:' + str(len(self.black_list[p])) + ']\n')
                    sys.exit(1)
                if len(input_variable[1]) == 1 and len(my_alt) == 1:
                    if self.black_list[p][my_var[0]]:
                        continue
                    self.snp_list[p].append(my_var)
                    self.black_list[p][my_var[0]] = 2
                else:
                    indel_failed = False
                    for k in range(my_var[0], my_var[0] + in_len):
                        if k >= len(self.black_list[p]):
                            indel_failed = True
                            continue
                        if self.black_list[p][k]:
                            indel_failed = True
                            continue
                    if indel_failed:
                        continue
                    for k in range(my_var[0], my_var[0] + in_len):
                        self.black_list[p][k] = 1
                    self.indel_list[p].append(my_var)

    def random_mutations(self):

        # add random indels
        all_indels = [[]] * len(self.sequences)
        for i in range(self.ploidy):
            for j in range(self.indels_to_add[i]):
                # insert homozygous indel
                if random.random() <= self.models[i][1]:
                    which_ploid = range(self.ploidy)
                # insert heterozygous indel
                else:
                    which_ploid = [self.ploid_mut_prior.sample()]

                # try to find suitable places to insert indels
                event_pos = -1
                for attempt in range(MAX_ATTEMPTS):
                    event_pos = random.randint(self.win_buffer, self.seq_len - 1)
                    for p in which_ploid:
                        if self.black_list[p][event_pos]:
                            event_pos = -1
                    if event_pos != -1:
                        break
                if event_pos == -1:
                    continue

                # insertion
                if random.random() <= self.models[i][3]:
                    in_len = self.models[i][4].sample()
                    # sequence content of random insertions is uniformly random (change this later, maybe)
                    in_seq = ''.join([random.choice(ALLOWED_NUCL) for _ in range(in_len)])
                    ref_nucl = self.sequences[i][event_pos]
                    my_indel = (event_pos, ref_nucl, ref_nucl + in_seq)
                # deletion
                else:
                    in_len = self.models[i][5].sample()
                    # skip if deletion too close to boundary
                    if event_pos + in_len + 1 >= len(self.sequences[i]):
                        continue
                    if in_len == 1:
                        in_seq = self.sequences[i][event_pos + 1]
                    else:
                        in_seq = str(self.sequences[i][event_pos + 1:event_pos + in_len + 1])
                    ref_nucl = self.sequences[i][event_pos]
                    my_indel = (event_pos, ref_nucl + in_seq, ref_nucl)

                # if event too close to boundary, skip. if event conflicts with other indel, skip.
                skip_event = False
                if event_pos + len(my_indel[1]) >= self.seq_len - self.win_buffer - 1:
                    skip_event = True
                if skip_event:
                    continue
                for p in which_ploid:
                    for k in range(event_pos, event_pos + in_len + 1):
                        if self.black_list[p][k]:
                            skip_event = True
                if skip_event:
                    continue

                for p in which_ploid:
                    for k in range(event_pos, event_pos + in_len + 1):
                        self.black_list[p][k] = 1
                    all_indels[p].append(my_indel)

        # add random snps
        all_snps = [[]] * self.ploidy
        for i in range(self.ploidy):
            for j in range(self.snps_to_add[i]):
                # insert homozygous SNP
                if random.random() <= self.models[i][1]:
                    which_ploid = range(self.ploidy)
                # insert heterozygous SNP
                else:
                    which_ploid = [self.ploid_mut_prior.sample()]

                # try to find suitable places to insert snps
                event_pos = -1
                for attempt in range(MAX_ATTEMPTS):
                    # based on the mutation model for the specified ploid, choose a SNP location based on trinuc bias
                    # (if there are multiple ploids, choose one at random)
                    if IGNORE_TRINUC:
                        event_pos = random.randint(self.win_buffer + 1, self.seq_len - 2)
                    else:
                        ploid_to_use = which_ploid[random.randint(0, len(which_ploid) - 1)]
                        event_pos = self.trinuc_bias[ploid_to_use].sample()
                    for p in which_ploid:
                        if self.black_list[p][event_pos]:
                            event_pos = -1
                    if event_pos != -1:
                        break
                if event_pos == -1:
                    continue

                ref_nucl = self.sequences[i][event_pos]
                context = str(self.sequences[i][event_pos - 1]) + str(self.sequences[i][event_pos + 1])
                # sample from tri-nucleotide substitution matrices to get SNP alt allele
                new_nucl = self.models[i][6][TRI_IND[context]][NUC_IND[ref_nucl]].sample()
                my_snp = (event_pos, ref_nucl, new_nucl)

                for p in which_ploid:
                    all_snps[p].append(my_snp)
                    self.black_list[p][my_snp[0]] = 2

        # combine random snps with inserted snps, remove any snps that overlap indels
        for p in range(len(all_snps)):
            all_snps[p].extend(self.snp_list[p])
            all_snps[p] = [n for n in all_snps[p] if self.black_list[p][n[0]] != 1]

        # MODIFY REFERENCE STRING: SNPS
        for i in range(len(all_snps)):
            temp = MutableSeq(self.sequences[i])
            for j in range(len(all_snps[i])):
                v_pos = all_snps[i][j][0]

                if all_snps[i][j][1] != temp[v_pos]:
                    print('\nError: Something went wrong!\n', all_snps[i][j], temp[v_pos], '\n')
                    print(all_snps[i][j])
                    print(self.sequences[i][v_pos])
                    sys.exit(1)
                else:
                    temp[v_pos] = all_snps[i][j][2]
            self.sequences[i] = Seq(temp)

        # organize the indels we want to insert
        for i in range(len(all_indels)):
            all_indels[i].extend(self.indel_list[i])
        all_indels_ins = [sorted([list(m) for m in n]) for n in all_indels]

        # MODIFY REFERENCE STRING: INDELS
        for i in range(len(all_indels_ins)):
            rolling_adj = 0
            # TODO looking to remove cigar string calculations
            temp_symbol_list = CigarString.string_to_list(str(len(self.sequences[i])) + "M")

            for j in range(len(all_indels_ins[i])):
                v_pos = all_indels_ins[i][j][0] + rolling_adj
                v_pos2 = v_pos + len(all_indels_ins[i][j][1])
                indel_length = len(all_indels_ins[i][j][2]) - len(all_indels_ins[i][j][1])
                rolling_adj += indel_length

                if all_indels_ins[i][j][1] != str(self.sequences[i][v_pos:v_pos2]):
                    print('\nError: Something went wrong!\n', all_indels_ins[i][j], [v_pos, v_pos2],
                          str(self.sequences[i][v_pos:v_pos2]), '\n')
                    sys.exit(1)
                else:
                    # alter reference sequence
                    self.sequences[i] = self.sequences[i][:v_pos] + Seq(all_indels_ins[i][j][2]) + \
                                        self.sequences[i][v_pos2:]
                    # notate indel positions for cigar computation
                    if indel_length > 0:
                        temp_symbol_list = temp_symbol_list[:v_pos + 1] + ['I'] * indel_length \
                                           + temp_symbol_list[v_pos2 + 1:]
                    elif indel_length < 0:
                        temp_symbol_list[v_pos + 1] = "D" * abs(indel_length) + "M"

            # pre-compute cigar strings
            for j in range(len(temp_symbol_list) - self.read_len):
                self.all_cigar[i].append(temp_symbol_list[j:j + self.read_len])

            # TODO do we need fm_pos at all?
            # create some data structures we will need later:
            # --- self.fm_pos[ploid][pos]: position of the left-most matching base (IN REFERENCE COORDINATES, i.e.
            #       corresponding to the unmodified reference genome)
            # --- self.fm_span[ploid][pos]: number of reference positions spanned by a read originating from
            #       this coordinate
            md_so_far = 0
            for j in range(len(temp_symbol_list)):
                self.fm_pos[i].append(md_so_far)
                # fix an edge case with deletions
                if 'D' in temp_symbol_list[j]:
                    self.fm_pos[i][-1] += temp_symbol_list[j].count('D')
                # compute number of ref matches for each read
                # This line gets hit a lot and is relatively slow. Might look for an improvement
                span_dif = len([n for n in temp_symbol_list[j: j + self.read_len] if 'M' in n])
                self.fm_span[i].append(self.fm_pos[i][-1] + span_dif)
                md_so_far += temp_symbol_list[j].count('M') + temp_symbol_list[j].count('D')

        # tally up all the variants we handled...
        count_dict = {}
        all_variants = [sorted(all_snps[i] + all_indels[i]) for i in range(self.ploidy)]
        for i in range(len(all_variants)):
            for j in range(len(all_variants[i])):
                all_variants[i][j] = tuple([all_variants[i][j][0] + self.x]) + all_variants[i][j][1:]
                t = tuple(all_variants[i][j])
                if t not in count_dict:
                    count_dict[t] = []
                count_dict[t].append(i)

        # TODO: combine multiple variants that happened to occur at same position into single vcf entry?

        output_variants = []
        for k in sorted(count_dict.keys()):
            output_variants.append(k + tuple([len(count_dict[k]) / float(self.ploidy)]))
            ploid_string = ['0' for _ in range(self.ploidy)]
            for k2 in [n for n in count_dict[k]]:
                ploid_string[k2] = '1'
            output_variants[-1] += tuple(['WP=' + '/'.join(ploid_string)])
        return output_variants

    def sample_read(self, sequencing_model, frag_len=None):

        # choose a ploid
        my_ploid = random.randint(0, self.ploidy - 1)

        # choose a random position within the ploid, and generate quality scores / sequencing errors
        reads_to_sample = []
        if frag_len is None:
            r_pos = self.coverage_distribution[my_ploid].sample()

            # sample read position and call function to compute quality scores / sequencing errors
            r_dat = self.sequences[my_ploid][r_pos:r_pos + self.read_len]
            (my_qual, my_errors) = sequencing_model.get_sequencing_errors(r_dat)
            reads_to_sample.append([r_pos, my_qual, my_errors, r_dat])

        else:
            r_pos1 = self.coverage_distribution[my_ploid][self.fraglen_ind_map[frag_len]].sample()

            # EXPERIMENTAL
            # coords_to_select_from = self.coverage_distribution[my_ploid][self.fraglens_ind_map[frag_len]].sample()
            # r_pos1 = random.randint(coords_to_select_from[0],coords_to_select_from[1])

            r_pos2 = r_pos1 + frag_len - self.read_len
            r_dat1 = self.sequences[my_ploid][r_pos1:r_pos1 + self.read_len]
            r_dat2 = self.sequences[my_ploid][r_pos2:r_pos2 + self.read_len]
            (my_qual1, my_errors1) = sequencing_model.get_sequencing_errors(r_dat1)
            (my_qual2, my_errors2) = sequencing_model.get_sequencing_errors(r_dat2, is_reverse_strand=True)
            reads_to_sample.append([r_pos1, my_qual1, my_errors1, r_dat1])
            reads_to_sample.append([r_pos2, my_qual2, my_errors2, r_dat2])

        # error format:
        # myError[i] = (type, len, pos, ref, alt)

        """
        examine sequencing errors to-be-inserted.
            - remove deletions that don't have enough bordering sequence content to "fill in"
            if error is valid, make the changes to the read data
        """
        read_out = []
        for read in reads_to_sample:
            try:
                my_cigar = self.all_cigar[my_ploid][read[0]]
            except IndexError:
                print('Index error when attempting to find cigar string.')
                print(my_ploid, len(self.all_cigar[my_ploid]), read[0])
                if frag_len is not None:
                    print((r_pos1, r_pos2))
                    print(frag_len, self.fraglen_ind_map[frag_len])
                sys.exit(1)
            total_d = sum([error[1] for error in read[2] if error[0] == 'D'])
            total_i = sum([error[1] for error in read[2] if error[0] == 'I'])
            avail_b = len(self.sequences[my_ploid]) - read[0] - self.read_len - 1

            # add buffer sequence to fill in positions that get deleted
            read[3] += self.sequences[my_ploid][read[0] + self.read_len:read[0] + self.read_len + total_d]
            # this is leftover code and a patch for a method that isn't used. There is probably a better
            # way to structure this than with a boolean
            first_time = True
            adj = 0
            sse_adj = [0 for _ in range(self.read_len + max(sequencing_model.err_p[3]))]
            any_indel_err = False

            # sort by letter (D > I > S) such that we introduce all indel errors before substitution errors
            # secondarily, sort by index
            arranged_errors = {'D': [], 'I': [], 'S': []}
            for error in read[2]:
                arranged_errors[error[0]].append((error[2], error))
            sorted_errors = []
            for k in sorted(arranged_errors.keys()):
                sorted_errors.extend([n[1] for n in sorted(arranged_errors[k])])

            skip_indels = False

            extra_cigar_val = []

            for error in sorted_errors:
                e_len = error[1]
                e_pos = error[2]
                if error[0] == 'D' or error[0] == 'I':
                    any_indel_err = True

                    if total_d > avail_b:  # if not enough bases to fill-in deletions, skip all indel erors
                        continue
                    if first_time:
                        # Again, this whole first time thing is a workaround for the previous
                        # code, which is simplified. May need to fix this all at some point
                        first_time = False
                        fill_to_go = total_d - total_i + 1
                        if fill_to_go > 0:
                            try:
                                extra_cigar_val = self.all_cigar[my_ploid][read[0] + fill_to_go][-fill_to_go:]
                            except IndexError:
                                # Applying the deletions we want requires going beyond region boundaries.
                                # Skip all indel errors
                                skip_indels = True

                    if skip_indels:
                        continue

                    # insert deletion error into read and update cigar string accordingly
                    if error[0] == 'D':
                        my_adj = sse_adj[e_pos]
                        pi = e_pos + my_adj
                        pf = e_pos + my_adj + e_len + 1
                        if str(read[3][pi:pf]) == str(error[3]):
                            read[3] = read[3][:pi + 1] + read[3][pf:]
                            my_cigar = my_cigar[:pi + 1] + my_cigar[pf:]
                            # weird edge case with del at very end of region. Make a guess and add a "M"
                            if pi + 1 == len(my_cigar):
                                my_cigar.append('M')

                            try:
                                my_cigar[pi + 1] = 'D' * e_len + my_cigar[pi + 1]
                            except IndexError:
                                print("Bug!! Index error on expanded cigar")
                                sys.exit(1)

                        else:
                            print('\nError, ref does not match alt while attempting to insert deletion error!\n')
                            sys.exit(1)
                        adj -= e_len
                        for i in range(e_pos, len(sse_adj)):
                            sse_adj[i] -= e_len

                    # insert insertion error into read and update cigar string accordingly
                    else:
                        my_adj = sse_adj[e_pos]
                        if str(read[3][e_pos + my_adj]) == error[3]:
                            read[3] = read[3][:e_pos + my_adj] + error[4] + read[3][e_pos + my_adj + 1:]
                            my_cigar = my_cigar[:e_pos + my_adj] + ['I'] * e_len + my_cigar[e_pos + my_adj:]
                        else:
                            print('\nError, ref does not match alt while attempting to insert insertion error!\n')
                            print('---', chr(read[3][e_pos + my_adj]), '!=', error[3])
                            sys.exit(1)
                        adj += e_len
                        for i in range(e_pos, len(sse_adj)):
                            sse_adj[i] += e_len

                else:  # substitution errors, much easier by comparison...
                    if str(read[3][e_pos + sse_adj[e_pos]]) == error[3]:
                        temp = MutableSeq(read[3])
                        temp[e_pos + sse_adj[e_pos]] = error[4]
                        read[3] = Seq(temp)
                    else:
                        print('\nError, ref does not match alt while attempting to insert substitution error!\n')
                        sys.exit(1)

            if any_indel_err:
                if len(my_cigar):
                    my_cigar = (my_cigar + extra_cigar_val)[:self.read_len]

                read[3] = read[3][:self.read_len]

            read_out.append([self.fm_pos[my_ploid][read[0]], my_cigar, read[3], str(read[1])])

        # read_out[i] = (pos, cigar, read_string, qual_string)
        return read_out
