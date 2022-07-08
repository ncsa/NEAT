import pickle, gzip
from neat.models import MutationModel
from neat.common import *


default_model = MutationModel(avg_mut_rate=0.001,
                              homozygous_freq=0.01,
                              insertion_chance=0.03,
                              deletion_chance=0.02,
                              substitution_chance=0.95,
                              trinuc_trans_matrices=trinuc_trans_matrices,
                              trinuc_trans_bias=DEFAULT_MUT_MODEL_TRINUC_BIAS,
                              insertion_lengths=DEFAULT_MUT_MODEL_INS_LENGTH_VALUES,
                              insertion_weights=DEFAULT_MUT_MODEL_INS_LENGTH_WEIGHTS,
                              deletion_lengths=DEFAULT_MUT_MODEL_DEL_LENGTH_VALUES,
                              deletion_weights=DEFAULT_MUT_MODEL_DEL_LENGTH_WEIGHTS,
                              is_cancer=False,)

default_cancer_model = MutationModel(avg_mut_rate=0.002,
                                     homozygous_freq=0.2,
                                     insertion_chance=0.03,
                                     deletion_chance=0.07,
                                     substitution_chance=0.9,
                                     is_cancer=False,
                                     trinuc_trans_matrices=trinuc_trans_matrices,
                                     trinuc_trans_bias=DEFAULT_MUT_MODEL_TRINUC_BIAS,
                                     insertion_lengths=DEFAULT_MUT_MODEL_INS_LENGTH_VALUES,
                                     insertion_weights=DEFAULT_MUT_MODEL_INS_LENGTH_WEIGHTS,
                                     deletion_lengths=DEFAULT_MUT_MODEL_DEL_LENGTH_VALUES,
                                     deletion_weights=DEFAULT_MUT_MODEL_DEL_LENGTH_WEIGHTS)

default_model_file = "/home/joshfactorial/NEAT/neat/models/defaults/default_mutation_model.pickle.gz"
with gzip.open(default_model_file, 'wt') as default:
    pickle.dump(default_model, default)

default_cancer_file = "/home/joshfactorial/NEAT/neat/models/defaults/default_cancer_mutation_model.pickle.gz"
with gzip.open(default_cancer_file, 'wt') as default_cancer:
    pickle.dump(default_cancer_model, default_cancer)

